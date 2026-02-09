from fastapi import APIRouter, Depends, HTTPException, Request
from fastapi.responses import RedirectResponse
from fastapi.templating import Jinja2Templates
from sqlalchemy.orm import Session
from .. import models, database, auth
from datetime import datetime

router = APIRouter(tags=["Dashboard"])
templates = Jinja2Templates(directory="app/templates")

# Mock Data (Moved to global scope for access by multiple routes)
MOCK_IND_REPORTS = [
    {
        "id": 101,
        "title": "IND Submission - Oncology Phase I",
        "date": datetime.now().strftime("%b %d, %Y"),
        "status": "Ready",
        "thumbnail_url": "https://via.placeholder.com/300x200/00487E/ffffff?text=Report+101",
        "description": "Comprehensive Investigational New Drug application for the novel oncology compound, including PK/PD modeling results and safety profile analysis.",
    },
    {
        "id": 102,
        "title": "Toxicity Profile Analysis",
        "date": datetime.now().strftime("%b %d, %Y"),
        "status": "In Progress",
        "thumbnail_url": "https://via.placeholder.com/300x200/00487E/ffffff?text=Report+102",
        "description": "Detailed analysis of potential toxicity risks based on in-silico predictions and early-stage in-vitro data.",
    },
    {
        "id": 103,
        "title": "FDA Safety Review Draft",
        "date": datetime.now().strftime("%b %d, %Y"),
        "status": "Draft",
        "thumbnail_url": "https://via.placeholder.com/300x200/00487E/ffffff?text=Report+103",
        "description": "Preliminary draft of the safety review section for FDA submission, adhering to the latest regulatory guidelines.",
    },
    {
        "id": 104,
        "title": "Clinical Trial Protocol v2",
        "date": datetime.now().strftime("%b %d, %Y"),
        "status": "Review",
        "thumbnail_url": "https://via.placeholder.com/300x200/00487E/ffffff?text=Report+104",
        "description": "Updated clinical trial protocol incorporating feedback from the initial review board and optimization results.",
    },
]


@router.get("/dashboard")
async def dashboard(request: Request, db: Session = Depends(database.get_db)):
    # Auth Check
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    # 1. PROJECTS Section (My Projects)
    projects = (
        db.query(models.Project)
        .filter(
            models.Project.owner_id == current_user.id,
            models.Project.is_deleted == False,
        )
        .all()
    )

    # 2. RESULTS Section (Grouped by Project)
    # Fetch projects that have at least one cohort (Simulation Result)
    # Using distinct to avoid duplicates if multiple cohorts exist
    results = (
        db.query(models.Project)
        .join(models.Cohort)
        .filter(
            models.Project.owner_id == current_user.id,
            models.Project.is_deleted == False,
        )
        .distinct()
        .all()
    )

    # 3. IND REPORT Section (From Database)
    ind_reports = (
        db.query(models.INDReport)
        .filter(models.INDReport.is_deleted == False)
        .order_by(models.INDReport.created_at.desc())
        .all()
    )

    return templates.TemplateResponse(
        "dashboard.html",
        {
            "request": request,
            "user": current_user,
            "projects": projects,
            "results": results,
            "ind_reports": ind_reports,
            "sidebar_mode": "default",  # Wide mode
        },
    )


@router.get("/dashboard/projects")
async def view_all_projects(request: Request, db: Session = Depends(database.get_db)):
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")
    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    projects = (
        db.query(models.Project)
        .filter(
            models.Project.owner_id == current_user.id,
            models.Project.is_deleted == False,
        )
        .all()
    )

    return templates.TemplateResponse(
        "projects_list.html",
        {
            "request": request,
            "user": current_user,
            "projects": projects,
            "sidebar_mode": "default",
        },
    )


@router.get("/dashboard/results")
async def view_all_results(request: Request, db: Session = Depends(database.get_db)):
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")
    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    results = (
        db.query(models.Project)
        .join(models.Cohort)
        .filter(
            models.Project.owner_id == current_user.id,
            models.Project.is_deleted == False,
        )
        .distinct()
        .all()
    )

    return templates.TemplateResponse(
        "results_list.html",
        {
            "request": request,
            "user": current_user,
            "results": results,
            "sidebar_mode": "default",
        },
    )


@router.get("/dashboard/ind-reports")
async def view_all_ind_reports(
    request: Request, db: Session = Depends(database.get_db)
):
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")
    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    ind_reports = (
        db.query(models.INDReport)
        .filter(models.INDReport.is_deleted == False)
        .order_by(models.INDReport.created_at.desc())
        .all()
    )
    # Need results for the modal if present
    results = (
        db.query(models.Project)
        .join(models.Cohort)
        .filter(models.Project.owner_id == current_user.id)
        .distinct()
        .all()
    )

    return templates.TemplateResponse(
        "ind_reports_list.html",
        {
            "request": request,
            "user": current_user,
            "ind_reports": ind_reports,
            "results": results,
            "sidebar_mode": "default",
        },
    )


@router.get("/dashboard/reports/{report_id}")
async def view_report_detail(
    request: Request, report_id: int, db: Session = Depends(database.get_db)
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

    # Find the report from Database
    report = db.query(models.INDReport).filter(models.INDReport.id == report_id).first()

    if not report:
        raise HTTPException(status_code=404, detail="Report not found")

    # Pass form_data from meta_data for pre-filling form fields
    form_data = report.meta_data if report.meta_data else {}

    # If report is completed, redirect to FDA 1571 PDF preview
    if report.status == "Completed" and form_data:
        from urllib.parse import urlencode
        params = {
            "project_id": report.project_id or "",
            "applicant_name": form_data.get("applicant_name", ""),
            "applicant_address": form_data.get("applicant_address", ""),
            "applicant_phone": form_data.get("applicant_phone", ""),
            "pi_name": form_data.get("pi_name", ""),
            "drug_name": form_data.get("drug_name", ""),
            "smiles": form_data.get("smiles", ""),
            "indication": form_data.get("indication", ""),
            "dose_mg": form_data.get("dose_mg", ""),
            "cmax": form_data.get("cmax", ""),
            "auc": form_data.get("auc", ""),
            "t_half": form_data.get("t_half", ""),
            "vss": form_data.get("vss", ""),
            "herg_margin": form_data.get("herg_margin", ""),
            "hepato_risk": form_data.get("hepato_risk", ""),
            "risk_category": form_data.get("risk_category", ""),
        }
        query_string = urlencode({k: v for k, v in params.items() if v})
        return RedirectResponse(url=f"/ind-generator/export-fda1571?{query_string}")

    return templates.TemplateResponse(
        "ind_report_detail.html",
        {
            "request": request,
            "user": current_user,
            "report": report,
            "form_data": form_data,
            "sidebar_mode": "compact",
        },
    )


@router.get("/dashboard/share")
async def share_summary(request: Request, db: Session = Depends(database.get_db)):
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    # Get shared projects (showing all for demo)
    shared_items = db.query(models.Project).all()

    return templates.TemplateResponse(
        "share_summary.html",
        {
            "request": request,
            "user": current_user,
            "shared_items": shared_items,
            "sidebar_mode": "default",
        },
    )


@router.get("/dashboard/users")
async def user_management(request: Request, db: Session = Depends(database.get_db)):
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    # Get all users
    users = db.query(models.User).all()

    return templates.TemplateResponse(
        "user_management.html",
        {
            "request": request,
            "user": current_user,
            "users": users,
            "sidebar_mode": "default",
        },
    )


@router.delete("/dashboard/reports/{report_id}")
async def delete_ind_report(
    request: Request, report_id: int, db: Session = Depends(database.get_db)
):
    token = request.cookies.get("access_token")
    if not token:
        raise HTTPException(status_code=401, detail="Not authenticated")

    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        raise HTTPException(status_code=401, detail="Not authenticated")

    report = db.query(models.INDReport).filter(models.INDReport.id == report_id).first()
    if not report:
        raise HTTPException(status_code=404, detail="Report not found")

    # Check ownership via project
    if report.project and report.project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Soft Delete
    report.is_deleted = True
    db.commit()

    return {"message": "Report moved to deleted items"}


@router.get("/dashboard/deleted")
async def view_deleted_objects(
    request: Request, db: Session = Depends(database.get_db)
):
    token = request.cookies.get("access_token")
    if not token:
        return RedirectResponse(url="/login")
    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        return RedirectResponse(url="/login")

    deleted_projects = (
        db.query(models.Project)
        .filter(
            models.Project.owner_id == current_user.id,
            models.Project.is_deleted == True,
        )
        .all()
    )

    deleted_reports = (
        db.query(models.INDReport)
        .join(models.Project)
        .filter(
            models.Project.owner_id == current_user.id,
            models.INDReport.is_deleted == True,
        )
        .all()
    )

    return templates.TemplateResponse(
        "deleted_objects.html",
        {
            "request": request,
            "user": current_user,
            "deleted_projects": deleted_projects,
            "deleted_reports": deleted_reports,
            "sidebar_mode": "default",
        },
    )


@router.post("/dashboard/restore/{type}/{id}")
async def restore_object(
    request: Request, type: str, id: int, db: Session = Depends(database.get_db)
):
    token = request.cookies.get("access_token")
    if not token:
        raise HTTPException(status_code=401)
    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        raise HTTPException(status_code=401)

    obj = None
    if type == "project":
        obj = (
            db.query(models.Project)
            .filter(models.Project.id == id, models.Project.owner_id == current_user.id)
            .first()
        )
    elif type == "report":
        obj = (
            db.query(models.INDReport)
            .join(models.Project)
            .filter(
                models.INDReport.id == id, models.Project.owner_id == current_user.id
            )
            .first()
        )

    if not obj:
        raise HTTPException(status_code=404, detail="Not found")

    obj.is_deleted = False
    db.commit()
    return {"status": "success"}


@router.delete("/dashboard/permanent/{type}/{id}")
async def permanent_delete_object(
    request: Request, type: str, id: int, db: Session = Depends(database.get_db)
):
    token = request.cookies.get("access_token")
    if not token:
        raise HTTPException(status_code=401)
    try:
        scheme, _, param = token.partition(" ")
        current_user = await auth.get_current_user(token=param, db=db)
    except:
        raise HTTPException(status_code=401)

    obj = None
    if type == "project":
        obj = (
            db.query(models.Project)
            .filter(models.Project.id == id, models.Project.owner_id == current_user.id)
            .first()
        )
    elif type == "report":
        obj = (
            db.query(models.INDReport)
            .join(models.Project)
            .filter(
                models.INDReport.id == id, models.Project.owner_id == current_user.id
            )
            .first()
        )

    if not obj:
        raise HTTPException(status_code=404, detail="Not found")

    db.delete(obj)
    db.commit()
    return {"status": "success"}
