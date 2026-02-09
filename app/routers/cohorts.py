from fastapi import APIRouter, Depends, HTTPException, Request, Form
from sqlalchemy.orm import Session
from fastapi.templating import Jinja2Templates
from fastapi.responses import RedirectResponse
import json

from .. import models, database, auth

router = APIRouter(
    prefix="/projects/{project_id}/cohorts",
    tags=["cohorts"],
    responses={404: {"description": "Not found"}},
)

templates = Jinja2Templates(directory="app/templates")


def get_db():
    db = database.SessionLocal()
    try:
        yield db
    finally:
        db.close()


@router.get("/new")
async def create_cohort_form(
    request: Request,
    project_id: int,
    pred_id: int = None,
    db: Session = Depends(database.get_db),
):
    # Auth Check (Cookie based)
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    project = db.query(models.Project).filter(models.Project.id == project_id).first()
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    # If pred_id is provided, fetch the prediction to pre-fill the form
    prediction = None
    if pred_id:
        prediction = (
            db.query(models.Prediction).filter(models.Prediction.id == pred_id).first()
        )

    # Prepare pre-fill data safely
    prefill_data = {"cl": "", "vd": "", "mw": "", "fu": "0.9"}
    if (
        prediction
        and prediction.result_json
        and isinstance(prediction.result_json, dict)
    ):
        pk = prediction.result_json.get("pk") or {}
        props = prediction.result_json.get("properties") or {}

        # CL: Try linear, then log
        cl_val = pk.get("human_CL_mL_min_kg_linear")
        if cl_val is None:
            cl_log = pk.get("human_CL_mL_min_kg")
            if cl_log is not None:
                try:
                    cl_val = 10 ** float(cl_log)
                except (ValueError, TypeError):
                    cl_val = None
        if cl_val is not None:
            prefill_data["cl"] = cl_val

        # VD: Try linear, then log
        vd_val = pk.get("human_VDss_L_kg_linear")
        if vd_val is None:
            vd_log = pk.get("human_VDss_L_kg")
            if vd_log is not None:
                try:
                    vd_val = 10 ** float(vd_log)
                except (ValueError, TypeError):
                    vd_val = None
        if vd_val is not None:
            prefill_data["vd"] = vd_val

        # MW
        if props.get("mw"):
            prefill_data["mw"] = props.get("mw")

        # FU
        if pk.get("human_fup"):
            prefill_data["fu"] = pk.get("human_fup")

    # Fetch existing cohorts for the list view (Defer heavy JSON results)
    from sqlalchemy.orm import defer

    cohorts = (
        db.query(models.Cohort)
        .filter(models.Cohort.project_id == project_id)
        .options(defer(models.Cohort.results_json))
        .order_by(models.Cohort.created_at.desc())
        .limit(20)
        .all()
    )

    return templates.TemplateResponse(
        "create_cohort.html",
        {
            "request": request,
            "project": project,
            "user": current_user,
            "prediction": prediction,
            "prefill_data": prefill_data,
            "cohorts": cohorts,
            "sidebar_mode": "compact",
            "is_compact": True,
        },
    )


from .. import inference


@router.post("")
async def create_cohort(
    request: Request,
    project_id: int,
    pop: str = Form(...),
    n_subjects: int = Form(...),
    age_min: int = Form(...),
    age_max: int = Form(...),
    female_ratio: float = Form(...),
    cl: float = Form(...),
    vd: float = Form(...),
    mw: float = Form(...),
    fu: float = Form(...),
    dose_mg: float = Form(100.0),
    phase: int = Form(1),  # Default to Phase 1
    # Advanced Settings (New)
    enzyme: str = Form("CYP2D6"),
    pheno_pm: float = Form(5.0),
    pheno_im: float = Form(25.0),
    pheno_em: float = Form(65.0),
    pheno_um: float = Form(5.0),
    bmi_min: float = Form(18.5),
    bmi_max: float = Form(30.0),
    bmi_dist: str = Form("normal"),
    renal_normal: float = Form(100.0),
    renal_mild: float = Form(0.0),
    renal_mod: float = Form(0.0),
    renal_sev: float = Form(0.0),
    hepatic_grade: str = Form("normal"),
    concomitant_meds: str = Form(""),
    db: Session = Depends(database.get_db),
):
    # Auth Check
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")
    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    project = db.query(models.Project).filter(models.Project.id == project_id).first()
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    # 1. Run Simulation
    # Note: In production, this should be a background task (Celery/Redis Queue)
    # For prototype, we run it synchronously (might take 5-10s)

    print(f"Starting Simulation for Project {project_id}...")
    results = inference.run_cohort_simulation(
        cl_ml_min_kg=cl,
        vd_l_kg=vd,
        mw_g_mol=mw,
        fup=fu,
        pop=pop,
        n_subjects=n_subjects,
        age_min=age_min,
        age_max=age_max,
        female_ratio=female_ratio,
        dose_mg=dose_mg,
    )

    if isinstance(results, dict) and "error" in results:
        # Handle error (flash message ideally, but simpler here)
        print(f"Simulation Failed: {results['error']}")
        # For now redirect back with error in query param?
        return RedirectResponse(
            url=f"/projects/{project_id}/cohorts/new?error=sim_failed", status_code=303
        )

    # 2. Save to DB
    # Generate a name
    count = (
        db.query(models.Cohort).filter(models.Cohort.project_id == project_id).count()
    )
    name = f"Cohort #{count + 1} ({pop}, n={n_subjects})"

    # Structure Advanced Settings
    advanced_settings = {
        "genetic": {
            "enzyme": enzyme,
            "distribution": {
                "PM": pheno_pm,
                "IM": pheno_im,
                "EM": pheno_em,
                "UM": pheno_um,
            },
        },
        "body": {"bmi_range": [bmi_min, bmi_max], "distribution": bmi_dist},
        "organ": {
            "renal": {
                "normal": renal_normal,
                "mild": renal_mild,
                "mod": renal_mod,
                "sev": renal_sev,
            },
            "hepatic": hepatic_grade,
        },
        "ddi": {
            # Parse comma-separated string if present
            "meds": [m.strip() for m in concomitant_meds.split(",") if m.strip()]
        },
    }

    demographics = {
        "pop": pop,
        "n": n_subjects,
        "age": [age_min, age_max],
        "female_ratio": female_ratio,
        "phase": phase,
        "advanced": advanced_settings,
    }
    drug_params = {"cl": cl, "vd": vd, "mw": mw, "fu": fu, "dose_mg": dose_mg}

    # We might want to compress results or just store full JSON (sqlite handles JSON)
    new_cohort = models.Cohort(
        project_id=project_id,
        name=name,
        description=f"Virtual cohort simulation for {pop}",
        demographics=demographics,
        drug_params=drug_params,
        results_json=results,
    )

    db.add(new_cohort)
    db.commit()
    db.refresh(new_cohort)

    # Redirect to the new cohort detail page (Phase specific views will be handled in Phase 2)
    return RedirectResponse(
        url=f"/projects/{project_id}/cohorts/{new_cohort.id}", status_code=303
    )


from ..inference import run_cohort_simulation, run_dose_response_analysis


@router.get("/{cohort_id}")
async def get_cohort_detail(
    request: Request,
    project_id: int,
    cohort_id: int,
    db: Session = Depends(database.get_db),
):
    # Auth Check
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")
    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    project = db.query(models.Project).filter(models.Project.id == project_id).first()
    cohort = db.query(models.Cohort).filter(models.Cohort.id == cohort_id).first()

    if not cohort:
        raise HTTPException(status_code=404, detail="Cohort not found")

    # Extract PK statistics from mPBPK simulation results
    results = cohort.results_json or []
    avg_cmax = 0.0
    avg_auc = 0.0
    pk_results = []
    
    # Safe extraction of drug_params
    drug_params = cohort.drug_params or {}
    mw = float(drug_params.get("mw", 150.0))

    if results and isinstance(results, list):
        completed_results = [
            r for r in results if isinstance(r, dict) and r.get("status") == "completed"
        ]
        if completed_results:
            print(
                f"DEBUG: First result keys: {completed_results[0].keys() if completed_results else 'None'}"
            )
            # Unit Conversion: nM -> µg/mL
            # Val(µg/mL) = Val(nM) * MW(g/mol) / 1,000,000
            cmax_values = [
                r.get("C_max", 0) * mw / 1000000.0 for r in completed_results
            ]
            # AUC is in nM*day. Convert to µg*hr/mL.
            # 1 day = 24 hr.
            auc_values = [
                r.get("AUC", 0) * 24 * mw / 1000000.0 for r in completed_results
            ]

            avg_cmax = sum(cmax_values) / len(cmax_values) if cmax_values else 0
            avg_auc = sum(auc_values) / len(auc_values) if auc_values else 0
            pk_results = completed_results
            print(
                f"DEBUG: avg_cmax={avg_cmax}, avg_auc={avg_auc}, n_completed={len(completed_results)}"
            )

    # Calculate Dose Response Data
    dose_response_data = []
    if cohort.drug_params:
        try:
            cl = cohort.drug_params.get("cl", 2.0)
            vd = cohort.drug_params.get("vd", 1.5)
            mw = cohort.drug_params.get("mw", 150.0)
            fu = cohort.drug_params.get("fu", 0.8)
            dose_mg = cohort.drug_params.get("dose_mg", 100.0)  # Default 100mg

            dose_response_data = run_dose_response_analysis(
                cl_ml_min_kg=float(cl),
                vd_l_kg=float(vd),
                mw_g_mol=float(mw),
                fup=float(fu),
                base_dose_mg=float(dose_mg),
            )
        except Exception as e:
            print(f"Dose Analysis Failed: {e}")

    # Get demographics for n_subjects
    demographics = cohort.demographics or {}
    n_subjects = demographics.get("n", demographics.get("n_subjects", 6))

    # Calculate Half-life (t½) from drug parameters
    # Formula: t½ = ln(2) × Vd / CL = 0.693 × Vd / CL
    # CL is in mL/min/kg, Vd is in L/kg
    # Convert CL: mL/min/kg → L/hr/kg: × 60 / 1000 = × 0.06
    avg_t12 = 0.0
    if cohort.drug_params:
        cl_ml_min_kg = cohort.drug_params.get("cl", 2.0)
        vd_l_kg = cohort.drug_params.get("vd", 1.5)
        cl_l_hr_kg = cl_ml_min_kg * 0.06  # Convert to L/hr/kg
        if cl_l_hr_kg > 0:
            avg_t12 = 0.693 * vd_l_kg / cl_l_hr_kg  # Result in hours

    return templates.TemplateResponse(
        "cohort_detail.html",
        {
            "request": request,
            "project": project,
            "cohort": cohort,
            "user": current_user,
            "dose_response_data": dose_response_data,
            "avg_cmax": round(avg_cmax, 1),
            "avg_auc": round(avg_auc, 1),
            "avg_t12": round(avg_t12, 1),
            "pk_results": pk_results,
            "n_subjects": n_subjects,
            "sidebar_mode": "compact",
            "is_compact": True,
        },
    )


@router.delete("/{cohort_id}")
async def delete_cohort(
    request: Request,
    project_id: int,
    cohort_id: int,
    db: Session = Depends(database.get_db),
):
    token = request.cookies.get("access_token")
    if not token:
        raise HTTPException(status_code=401, detail="Not authenticated")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        raise HTTPException(status_code=401, detail="Not authenticated")

    cohort = (
        db.query(models.Cohort)
        .filter(models.Cohort.id == cohort_id, models.Cohort.project_id == project_id)
        .first()
    )

    if not cohort:
        raise HTTPException(status_code=404, detail="Cohort not found")

    # Ensure the project belongs to the user
    if cohort.project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    db.delete(cohort)
    db.commit()

    return {"status": "success"}
