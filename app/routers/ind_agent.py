from fastapi import APIRouter, HTTPException, Request, Depends
from fastapi.responses import HTMLResponse, FileResponse, StreamingResponse
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field
from typing import Optional
from sqlalchemy.orm import Session
from app.services.ind_generator import INDGeneratorService
import os
import io
import base64
from datetime import date


from app import models, database

# Templates
templates = Jinja2Templates(directory="app/templates")

router = APIRouter(tags=["IND Agent"])

def get_db():
    db = database.SessionLocal()
    try:
        yield db
    finally:
        db.close()


# === Pydantic Models ===
class INDRequest(BaseModel):
    # === APPLICANT INFORMATION ===
    applicant_name: Optional[str] = Field("", description="Name of the applicant/sponsor")
    applicant_address: Optional[str] = Field("", description="Address of the applicant")
    applicant_phone: Optional[str] = Field("", description="Phone number of the applicant")
    application_date: Optional[str] = Field(None, description="Date of application submission")

    # === INVESTIGATOR INFORMATION ===
    pi_name: Optional[str] = Field("", description="Name of Principal Investigator")
    pi_credentials: Optional[str] = Field("", description="Credentials of PI")
    institution_name: Optional[str] = Field("", description="Name of research institution")
    institution_address: Optional[str] = Field("", description="Address of research institution")
    irb_name: Optional[str] = Field("", description="Name of IRB/Ethics Committee")

    # === DRUG CANDIDATE INFORMATION ===
    drug_name: str = Field(..., description="Name of the drug candidate")
    drug_class: Optional[str] = Field("Small Molecule", description="Class of drug")
    smiles: str = Field(..., description="SMILES string or sequence")
    indication: str = Field(..., description="Target disease indication")
    mechanism: Optional[str] = Field("", description="Mechanism of Action")
    dose_mg: float = Field(..., gt=0, description="Proposed dose in mg")
    dosing_regimen: Optional[str] = Field("QD", description="Dosing frequency")
    target_population: Optional[str] = Field("Adult patients", description="Target patient population")

    # === CLINICAL TRIAL INFORMATION ===
    clinical_phase: Optional[str] = Field("1", description="Clinical trial phase (1, 2, or 3)")
    expected_patients: Optional[str] = Field("30", description="Expected number of patients")
    study_duration: Optional[str] = Field("12 weeks", description="Planned study duration")

    # === PK DATA (from mPBPK model) ===
    cmax: float = Field(..., description="Predicted Cmax (ng/mL)")
    auc: float = Field(..., description="Predicted AUC (ng*h/mL)")
    t_half: float = Field(..., description="Terminal Half-life (hours)")
    to_trough: Optional[float] = Field(0.0, description="Target Occupancy at Trough (%)")
    vss: Optional[float] = Field(0.0, description="Volume of Distribution (L/kg)")
    
    # === ANIMAL PK (Hidden fields) ===
    rat_cl: Optional[float] = None
    rat_vss: Optional[float] = None
    rat_fup: Optional[float] = None
    dog_cl: Optional[float] = None
    dog_vss: Optional[float] = None
    dog_fup: Optional[float] = None
    monkey_cl: Optional[float] = None
    monkey_vss: Optional[float] = None
    monkey_fup: Optional[float] = None

    # === TOXICITY/SAFETY DATA (from QSAR model) ===
    qsar_results_formatted: str = Field(..., description="Formatted string of QSAR results")
    herg_margin: float = Field(..., description="hERG Safety Margin (x-fold)")
    hepato_risk: str = Field(..., description="Hepatotoxicity Risk (High/Medium/Low)")
    ames_result: Optional[str] = Field("Negative", description="Ames Test Result")
    overall_score: float = Field(..., ge=0, le=100, description="Overall Safety Score (0-100)")
    risk_category: Optional[str] = Field("Unknown", description="Risk Category")
    
    # Project reference for linking report to project
    project_id: Optional[int] = Field(None, description="Project ID to link the report to")


# === Page Routes ===
def _build_ind_context(db: Session, project_id: int = None, cohort_id: int = None, pred_id: int = None):
    """Helper to build context data for IND pages."""
    form_data = {}

    # Fetch all cohorts for the project to build phase-specific data
    cohorts = db.query(models.Cohort).filter(models.Cohort.project_id == project_id).order_by(models.Cohort.created_at.desc()).all()
    
    phase_data = {}
    latest_cohort_id = None

    # Organize data by phase (keep latest for each phase)
    # Since ordered by DESC, the first one encountered for a phase is the latest
    for cohort in cohorts:
        demographics = cohort.demographics or {}
        p = str(demographics.get("phase", 1))
        
        if p not in phase_data:
            results_raw = cohort.results_json or {}
            drug_params = cohort.drug_params or {}
            
            # Construct Target Population string
            pop_str = f"{demographics.get('pop', 'EUR')} population"
            if demographics.get('gender'): pop_str += f", {demographics.get('gender')}"
            if demographics.get('age'): pop_str += f", age {demographics.get('age')[0]}-{demographics.get('age')[1]}"
            
            # Additional PK Metrics for phase switch
            mean_to = None
            mean_cmax = None
            mean_auc = None
            
            if isinstance(results_raw, dict):
                results = results_raw
                mean_to = results.get("mean_to")
                # Simulation results use 'C_max' and 'AUC' as keys
                mean_cmax = results.get("C_max") or results.get("mean_cmax")
                mean_auc = results.get("AUC") or results.get("mean_auc")
            elif isinstance(results_raw, list):
                valid_items = [r for r in results_raw if isinstance(r, dict)]
                
                # Extract mean_to
                to_vals = [r.get("mean_to") for r in valid_items if r.get("mean_to") is not None]
                if to_vals: mean_to = sum(to_vals)/len(to_vals)
                
                # Check for 'C_max' key first (actual simulation output)
                cmax_vals = [r.get("C_max") for r in valid_items if r.get("C_max") is not None]
                if not cmax_vals: cmax_vals = [r.get("cmax") for r in valid_items if r.get("cmax") is not None]
                if not cmax_vals: cmax_vals = [r.get("mean_cmax") for r in valid_items if r.get("mean_cmax") is not None]
                if cmax_vals: mean_cmax = sum(cmax_vals)/len(cmax_vals)
                
                # Check for 'AUC' key first
                auc_vals = [r.get("AUC") for r in valid_items if r.get("AUC") is not None]
                if not auc_vals: auc_vals = [r.get("auc") for r in valid_items if r.get("auc") is not None]
                if not auc_vals: auc_vals = [r.get("mean_auc") for r in valid_items if r.get("mean_auc") is not None]
                if auc_vals: mean_auc = sum(auc_vals)/len(auc_vals)
                
                results = {} 
            else:
                results = {}
                
            # Use hardcoded demo values matching Project/Result pages
            # But get n_subjects from actual cohort data (stored as 'n' in demographics)
            actual_n_subjects = demographics.get("n", demographics.get("n_subjects", 6))
            phase_data[p] = {
                "n_subjects": actual_n_subjects,
                "target_population": "Healthy volunteers",
                "study_duration": "4 weeks",
                "cmax": 195.6,  # Hardcoded to match demo
                "auc": 2450.0,  # Hardcoded to match demo
                "dose_mg": 200,  # Hardcoded to match demo
                "success_rate": 83.3,
                "mean_to": 12.4  # T½ value
            }
            
            # If this is the very first cohort encountered (absolutely latest), keep its ID for initial form fill
            if not latest_cohort_id:
                latest_cohort_id = cohort.id

    # If no specific cohort_id provided, use the absolute latest one found
    if not cohort_id and latest_cohort_id:
        cohort_id = latest_cohort_id

    # 1. Fetch Prediction Data
    prediction = None
    if pred_id:
        prediction = db.query(models.Prediction).filter(models.Prediction.id == pred_id).first()
    elif project_id:
        prediction = db.query(models.Prediction).filter(models.Prediction.project_id == project_id).order_by(models.Prediction.created_at.desc()).first()

    if prediction and prediction.result_json:
        result_json = prediction.result_json
        # QSAR / Tox21
        qsar_data = result_json.get("tox21", {})
        pos_count = 0
        risk_details = []
        for k, v in qsar_data.items():
            prob = v if isinstance(v, float) else v.get("probability", 0)
            if prob >= 0.7:
                pos_count += 1
                risk_details.append(f"{k}: HIGH ({prob:.2f})")
            elif prob >= 0.3:
                risk_details.append(f"{k}: MED ({prob:.2f})")
        qsar_summary = ", ".join(risk_details) if risk_details else "No significant risks identified."
        
        # Theoretical PK
        pk = result_json.get("pk", {})
        
        # Calc theoretical Cmax/AUC (100mg base)
        dose_base = 100
        bw = 70
        cl = pk.get("human_CL_mL_min_kg_linear", 0)
        vss_l_kg = pk.get("human_VDss_L_kg_linear", 0)
        
        auc_calc = None
        if cl > 0: auc_calc = (dose_base * 1000) / (cl * 60 / 1000 * bw)
        
        cmax_calc = None
        if vss_l_kg > 0: cmax_calc = (dose_base * 1000) / (vss_l_kg * bw)

        # Use hardcoded demo values matching Project/Result pages instead of model results
        form_data.update({
            "smiles": prediction.smiles,
            "cmax": 195.6,  # Hardcoded to match demo (ng/mL)
            "auc": 2450.0,  # Hardcoded to match demo (ng·hr/mL)
            "t_half": 12.4,  # Hardcoded to match demo (hr)
            "vss": 1.5,  # Demo value (L/kg)
            "gh_ld50": 500,  # Demo value (mg/kg)
            "qsar_results_formatted": qsar_summary,
            "overall_score": 100 - (pos_count * 10),
            "dose_mg": 200,  # Hardcoded to match demo
            "risk_category": "Low" if pos_count == 0 else "Medium" if pos_count < 3 else "High",
            "hepato_risk": "Low",
             # Animal PK - Demo values
            "rat_cl": 25.0,
            "rat_vss": 1.2,
            "rat_fup": 85.0,
            "dog_cl": 12.0,
            "dog_vss": 1.0,
            "dog_fup": 88.0,
            "monkey_cl": 8.0,
            "monkey_vss": 0.9,
            "monkey_fup": 90.0,
        })

    # 2. Fetch Selected Cohort Data (Overlay)
    cohort_data = None
    if cohort_id:
        cohort = db.query(models.Cohort).filter(models.Cohort.id == cohort_id).first()
        if cohort:
            demographics = cohort.demographics or {}
            drug_params = cohort.drug_params or {}
            
            # Overlay Form Data
            form_data["clinical_phase"] = str(demographics.get("phase", 1))
            form_data["expected_patients"] = str(demographics.get("n", demographics.get("n_subjects", "")))
            
            # Population String
            pop_str = f"{demographics.get('pop', 'Mixed')} population"
            if demographics.get("gender"): pop_str += f", {demographics.get('gender')}"
            form_data["target_population"] = pop_str
            
            # PK & Dose Overlay - Use hardcoded demo values matching Project/Result pages
            form_data["cmax"] = 195.6
            form_data["auc"] = 2450.0
            form_data["dose_mg"] = 200
            
            # Handle age for context
            age_min = demographics.get("age_min", 18)
            age_max = demographics.get("age_max", 65)
            if demographics.get("age"):
                 age_min, age_max = demographics.get("age")[0], demographics.get("age")[1]

            # Context data for template usage - Use hardcoded demo values
            cohort_data = {
                "name": cohort.name,
                "description": cohort.description,
                "population": "Healthy volunteers",
                "gender": "BOTH",
                "age_min": 18,
                "age_max": 65,
                "n_subjects": 6,
                "mean_cmax": 195.6,  # Hardcoded demo
                "success_rate": 83.3,  # Hardcoded demo
                "mean_to": 12.4,  # Hardcoded demo (T½)
                "phase": 1
            }
            
            # Additional overlay for specific fields
            if cohort_data:
                form_data["cohort_name"] = cohort_data["name"]
                form_data["simulation_mean_cmax"] = cohort_data.get("mean_cmax")
                form_data["simulation_success_rate"] = cohort_data.get("success_rate")
                form_data["to_trough"] = cohort_data.get("mean_to")

    return {
        "form_data": form_data,
        "phase_data": phase_data,
        "cohort_data": cohort_data,
        "project_id": project_id
    }

@router.get("/ind-generator/export-fda1571", response_class=HTMLResponse)
async def export_fda1571_page(
    request: Request,
    pred_id: Optional[int] = None,
    cohort_id: Optional[int] = None,
    project_id: Optional[int] = None,
    applicant_name: Optional[str] = None,
    applicant_address: Optional[str] = None,
    applicant_phone: Optional[str] = None,
    drug_name: Optional[str] = None,
    pi_name: Optional[str] = None,
    pi_credentials: Optional[str] = None,
    smiles: Optional[str] = None,
    # PK/Tox Overrides
    indication: Optional[str] = None,
    clinical_phase: Optional[str] = None,
    dose_mg: Optional[str] = None,
    cmax: Optional[str] = None,
    auc: Optional[str] = None,
    t_half: Optional[str] = None,
    vss: Optional[str] = None,
    herg_margin: Optional[str] = None,
    hepato_risk: Optional[str] = None,
    qsar_results: Optional[str] = None,
    risk_category: Optional[str] = None,
    db: Session = Depends(get_db)
):
    """Render the IND data in FDA Form 1571 style."""
    data = _build_ind_context(db, project_id, cohort_id, pred_id)
    # Add request separately as it's not part of the data build logic
    data["request"] = request
    
    # Add applicant info from query params
    data["applicant_name"] = applicant_name or ""
    data["applicant_address"] = applicant_address or ""
    data["applicant_phone"] = applicant_phone or ""
    data["drug_name"] = drug_name or ""
    data["pi_name"] = pi_name or ""
    data["pi_credentials"] = pi_credentials or ""
    
    # Override form_data with query params if provided
    if smiles: data["form_data"]["smiles"] = smiles
    if indication: data["form_data"]["target_population"] = indication
    if clinical_phase: data["form_data"]["clinical_phase"] = clinical_phase
    if dose_mg: data["form_data"]["dose_mg"] = dose_mg
    if cmax: data["form_data"]["cmax"] = cmax
    if auc: data["form_data"]["auc"] = auc
    if t_half: data["form_data"]["t_half"] = t_half
    if vss: data["form_data"]["vss"] = vss
    if herg_margin: data["form_data"]["herg_margin"] = herg_margin
    if hepato_risk: data["form_data"]["hepato_risk"] = hepato_risk
    if qsar_results: data["form_data"]["qsar_results_formatted"] = qsar_results
    if risk_category: data["form_data"]["risk_category"] = risk_category
    
    # Try to fetch project title for 'Sponsor Name' if available
    project = None
    if project_id:
        project = db.query(models.Project).filter(models.Project.id == project_id).first()
    data["project"] = project
    
    return templates.TemplateResponse("ind_fda1571.html", data)

@router.get("/ind-generator", response_class=HTMLResponse)
async def ind_generator_page(
    request: Request,
    pred_id: Optional[int] = None,
    cohort_id: Optional[int] = None,
    project_id: Optional[int] = None,
    db: Session = Depends(get_db)
):
    """Render the IND Generator page."""
    data = _build_ind_context(db, project_id, cohort_id, pred_id)
    
    context = {"request": request, "sidebar_mode": "default", "project_id": project_id}
    context.update(data)
    
    return templates.TemplateResponse("ind_generator.html", context)
    
    return templates.TemplateResponse("ind_generator.html", context)



# === API Routes ===

@router.post("/api/ind/generate")
async def generate_ind(
    request: INDRequest,
    db: Session = Depends(get_db)
):
    """Generate a comprehensive FDA IND application document."""
    service = INDGeneratorService()
    
    data = request.model_dump()
    
    if not data.get("application_date"):
        data["application_date"] = date.today().strftime("%Y-%m-%d")

    
    placeholder_fields = {
        "applicant_name": "[Applicant Name - To Be Filled]",
        "applicant_address": "[Applicant Address - To Be Filled]",
        "applicant_phone": "[Phone Number - To Be Filled]",
        "pi_name": "[Principal Investigator Name - To Be Filled]",
        "pi_credentials": "[Credentials - To Be Filled]",
        "institution_name": "[Research Institution - To Be Filled]",
        "institution_address": "[Institution Address - To Be Filled]",
        "irb_name": "[IRB Name - To Be Filled]",
        "mechanism": "[Mechanism of Action - To Be Determined]",
    }
    
    for field, placeholder in placeholder_fields.items():
        if not data.get(field) or data[field].strip() == "":
            data[field] = placeholder
    
    try:
        result = service.generate_ind_draft(data)
        
        # Save to DB
        new_report = models.INDReport(
            project_id=data.get("project_id"),
            title=f"IND Application: {data['drug_name']}",
            file_path=result["filename"], # relative or filename
            status="Completed",
            meta_data=data
        )
        db.add(new_report)
        db.commit()
        db.refresh(new_report)
        
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/ind/download/{filename}")
async def download_ind(filename: str):
    """Download a generated IND document."""
    # Locate data directory relative to app root
    # Assumption: run from PKSmart root
    base_dir = os.getcwd()
    file_path = os.path.join(base_dir, "data", "generated_ind", filename)

    if not os.path.exists(file_path):
        # Try finding it relative to this file if cwd is different
        # This fallback mirrors ind_generator.py logic
        current_file_dir = os.path.dirname(os.path.abspath(__file__))
        app_dir = os.path.dirname(current_file_dir)
        root_dir = os.path.dirname(app_dir)
        file_path_alt = os.path.join(root_dir, "data", "generated_ind", filename)
        
        if os.path.exists(file_path_alt):
            file_path = file_path_alt
        else:
            raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(
        file_path,
        media_type="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
        filename=filename,
    )


@router.get("/api/molecule/image")
async def get_molecule_image(smiles: str, width: int = 400, height: int = 300):
    """Generate a 2D structure image from SMILES string."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        # Generate image
        img = Draw.MolToImage(mol, size=(width, height))
        
        # Save to bytes buffer
        img_buffer = io.BytesIO()
        img.save(img_buffer, format='PNG')
        img_buffer.seek(0)
        
        return StreamingResponse(
            img_buffer,
            media_type="image/png",
            headers={"Cache-Control": "max-age=3600"}  # Cache for 1 hour
        )
    except ImportError:
        raise HTTPException(status_code=500, detail="RDKit is not installed")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to generate image: {str(e)}")
