from fastapi import APIRouter, Depends, HTTPException, status, Request, Form
from fastapi.responses import RedirectResponse
from fastapi.templating import Jinja2Templates
from sqlalchemy.orm import Session
from typing import Optional
from .. import schemas, models, database, auth

router = APIRouter(tags=["Projects"])
templates = Jinja2Templates(directory="app/templates")


@router.get("/projects/{project_id}/results")
async def view_project_results(
    request: Request, project_id: int, db: Session = Depends(database.get_db)
):
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    project = (
        db.query(models.Project)
        .filter(
            models.Project.id == project_id, models.Project.owner_id == current_user.id
        )
        .first()
    )
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    # Get all cohorts for this project
    cohorts = (
        db.query(models.Cohort).filter(models.Cohort.project_id == project_id).all()
    )

    # Calculate Phase 1 Stats (Latest Cohort)
    p1_stats = {
        "dose": 0,
        "cmax": 0,
        "tmax": 0,
        "t12": 0,
        "auc": 0,
        "dlt_rate": 0,
        "n_subjects": 0,
    }

    # Data for PK Profile Charts
    pk_time_profile = []
    dose_comparison_data = []

    # Safe phase extraction helper
    def get_phase(c):
        demo = c.demographics or {}
        return str(demo.get("phase", "1"))

    p1_cohorts = [c for c in cohorts if get_phase(c) == "1"]

    # ... (Phase 2/3 existence check logic) ...
    # Check for existence of Phase 2 and Phase 3
    # Enable Phase 2/3 tabs by default to show demo data
    p2_cohorts = [c for c in cohorts if get_phase(c) == "2"]
    p3_cohorts = [c for c in cohorts if get_phase(c) == "3"]
    has_phase2 = len(p2_cohorts) > 0
    has_phase3 = len(p3_cohorts) > 0

    if p1_cohorts:
        # Get latest Phase 1 cohort
        latest_p1 = sorted(p1_cohorts, key=lambda x: x.id, reverse=True)[0]

        # Get Params (Safe Extraction)
        drug_params = latest_p1.drug_params or {}
        demographics = latest_p1.demographics or {}
        dose = float(drug_params.get("dose_mg") or 100.0)
        mw = float(drug_params.get("mw") or 150.0)
        cl = float(drug_params.get("cl") or 2.0)
        vd = float(drug_params.get("vd") or 1.5)
        fu = float(drug_params.get("fu") or 0.8)
        n_subjects = demographics.get("n", 6)

        p1_stats["dose"] = dose
        p1_stats["n_subjects"] = n_subjects

        # Calculate t1/2 (Same as before)
        cl_l_hr_kg = cl * 0.06
        if cl_l_hr_kg > 0:
            t12 = 0.693 * vd / cl_l_hr_kg
            p1_stats["t12"] = round(t12, 1)

        # Calculate Cmax, AUC from results
        results = latest_p1.results_json
        if results and isinstance(results, list):
            valid_res = [
                r
                for r in results
                if isinstance(r, dict) and r.get("status") == "completed"
            ]
            if valid_res:
                # 1. Update Stats
                cmax_vals = [
                    (r.get("C_max") or 0) * mw / 1000000.0 for r in valid_res
                ]  # µg/mL
                auc_vals = [
                    (r.get("AUC") or 0) * 24 * mw / 1000000.0 for r in valid_res
                ]  # µg*hr/mL

                if cmax_vals:
                    p1_stats["cmax"] = round(sum(cmax_vals) / len(cmax_vals), 1)
                if auc_vals:
                    p1_stats["auc"] = round(sum(auc_vals) / len(auc_vals), 1)

                # 2. Extract PK Time Profile (Spaghetti Plot)
                # Sample up to 10 subjects to keep chart readable and light
                import numpy as np

                for i, r in enumerate(valid_res[:10]):
                    if "time" in r and "C_plasma" in r:
                        # Downsample points if too many (e.g., take every 10th point)
                        times = np.array(r["time"]) * 24.0  # Days to Hours
                        concs = np.array(r["C_plasma"]) * mw  # nM to ng/mL

                        # Filter for reasonable number of points (e.g., 50 points)
                        indices = np.linspace(0, len(times) - 1, 50, dtype=int)

                        pk_time_profile.append(
                            {
                                "id": i + 1,
                                "time": times[indices].tolist(),
                                "conc": concs[indices].tolist(),
                            }
                        )

        # 3. Generate Dose Comparison Data (Simulate 50, 100, 150, RealDose)
        # Define doses: 50, 100, 150 + Actual Dose
        sample_doses = [50.0, 100.0, 150.0]
        if dose not in sample_doses:
            sample_doses.append(dose)
        sample_doses.sort()

        # Run Simulation
        from ..inference import run_dose_response_analysis

        dose_comparison_data = run_dose_response_analysis(
            cl_ml_min_kg=cl, vd_l_kg=vd, mw_g_mol=mw, fup=fu, compare_doses=sample_doses
        )

    # Calculate Phase 2 Stats
    p2_stats = {"n_subjects": 0, "dose": 0}
    p2_cohorts = [c for c in cohorts if get_phase(c) == "2"]
    if p2_cohorts:
        latest_p2 = sorted(p2_cohorts, key=lambda x: x.id, reverse=True)[0]
        demo_p2 = latest_p2.demographics or {}
        drug_p2 = latest_p2.drug_params or {}
        p2_stats["n_subjects"] = demo_p2.get("n", 0)
        p2_stats["dose"] = float(drug_p2.get("dose_mg") or 150.0)

    # Calculate Phase 3 Stats
    p3_stats = {"n_subjects": 0}
    p3_cohorts = [c for c in cohorts if get_phase(c) == "3"]
    if p3_cohorts:
        latest_p3 = sorted(p3_cohorts, key=lambda x: x.id, reverse=True)[0]
        demo_p3 = latest_p3.demographics or {}
        p3_stats["n_subjects"] = demo_p3.get("n", 0)

    # Check for existence of Phase 2 and Phase 3
    has_phase2 = len(p2_cohorts) > 0
    has_phase3 = len(p3_cohorts) > 0

    if p1_cohorts:
        # Get latest Phase 1 cohort
        latest_p1 = sorted(p1_cohorts, key=lambda x: x.id, reverse=True)[0]

        # Get Params (Safe Extraction)
        drug_params = latest_p1.drug_params or {}
        demographics = latest_p1.demographics or {}
        dose = drug_params.get("dose_mg") or 100.0
        mw = float(drug_params.get("mw") or 150.0)
        cl = float(drug_params.get("cl") or 2.0)
        vd = float(drug_params.get("vd") or 1.5)
        n_subjects = demographics.get("n", 6)

        p1_stats["dose"] = dose
        p1_stats["n_subjects"] = n_subjects

        # Calculate t1/2
        # CL (mL/min/kg) -> L/hr/kg: * 60 / 1000 = * 0.06
        cl_l_hr_kg = cl * 0.06
        if cl_l_hr_kg > 0:
            t12 = 0.693 * vd / cl_l_hr_kg
            p1_stats["t12"] = round(t12, 1)

        # Calculate Cmax, AUC from results
        results = latest_p1.results_json
        if results and isinstance(results, list):
            valid_res = [
                r
                for r in results
                if isinstance(r, dict) and r.get("status") == "completed"
            ]
            if valid_res:
                # Cmax (nM) -> µg/mL: Val(nM) * MW / 1,000,000
                cmax_vals = [(r.get("C_max") or 0) * mw / 1000000.0 for r in valid_res]
                # AUC (nM*day) -> µg*hr/mL: Val * 24 * MW / 1,000,000
                auc_vals = [
                    (r.get("AUC") or 0) * 24 * mw / 1000000.0 for r in valid_res
                ]

                if cmax_vals:
                    p1_stats["cmax"] = round(sum(cmax_vals) / len(cmax_vals), 1)
                if auc_vals:
                    p1_stats["auc"] = round(sum(auc_vals) / len(auc_vals), 1)

    # Calculate Phase 3 Stats
    p3_stats = {"n_subjects": 0}
    if p3_cohorts:
        latest_p3 = sorted(p3_cohorts, key=lambda x: x.id, reverse=True)[0]
        demo_p3 = latest_p3.demographics or {}
        p3_stats["n_subjects"] = demo_p3.get("n", 0)

    return templates.TemplateResponse(
        "project_results.html",
        {
            "request": request,
            "user": current_user,
            "project": project,
            "cohorts": cohorts,
            "p1_stats": p1_stats,
            "p2_stats": p2_stats,
            "p3_stats": p3_stats,
            "has_phase2": has_phase2,
            "has_phase3": has_phase3,
            "pk_time_profile": pk_time_profile,
            "dose_comparison_data": dose_comparison_data,
            "sidebar_mode": "compact",  # Compact sidebar as requested
            "sidebar_mode": "compact",  # Compact sidebar as requested
            "is_compact": True,
        },
    )


@router.post("/projects")
async def create_project(
    request: Request,
    title: str = Form(...),
    description: str = Form(None),
    db: Session = Depends(database.get_db),
):
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    new_project = models.Project(
        title=title, description=description, owner_id=current_user.id
    )
    db.add(new_project)
    db.commit()

    return RedirectResponse(url="/dashboard", status_code=303)


@router.get("/projects/{project_id}")
async def view_project(
    request: Request, project_id: int, db: Session = Depends(database.get_db)
):
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    project = (
        db.query(models.Project)
        .filter(
            models.Project.id == project_id, models.Project.owner_id == current_user.id
        )
        .first()
    )
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    return templates.TemplateResponse(
        "project_detail.html",
        {
            "request": request,
            "user": current_user,
            "project": project,
            "sidebar_mode": "compact",
            "is_compact": True,
        },
    )


@router.delete("/projects/{project_id}")
async def delete_project(
    request: Request, project_id: int, db: Session = Depends(database.get_db)
):
    token = request.cookies.get("access_token")
    if not token:
        raise HTTPException(status_code=401, detail="Not authenticated")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        raise HTTPException(status_code=401, detail="Not authenticated")

    # Find project and ensure ownership
    project = (
        db.query(models.Project)
        .filter(
            models.Project.id == project_id, models.Project.owner_id == current_user.id
        )
        .first()
    )

    if not project:
        raise HTTPException(
            status_code=404, detail="Project not found or permission denied"
        )

    # Soft Delete
    project.is_deleted = True
    db.commit()

    return {"message": "Project moved to deleted items"}
