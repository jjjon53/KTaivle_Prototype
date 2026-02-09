# 🏥 PKSmart: AI-Powered Clinical Trial Design Platform

Pyxidis의 AI 기반 임상시험 설계 및 독성 예측 플랫폼입니다.

---

## ✅ 주요 기능 (Implemented Features)

### 1. 🧪 QSAR 독성 예측 (Tox21)
- 12개 Tox21 엔드포인트에 대한 독성 확률 예측
- SMILES 입력으로 분자 구조 기반 분석
- 실시간 예측 결과 시각화

### 2. 💊 mPBPK 시뮬레이션
- **가상 코호트 생성**: 인구집단별 (EUR, EAS, AFR 등) 가상 피험자 생성
- **CYP2D6 대사 다형성 반영**: PharmVar 데이터 기반 대사 표현형 시뮬레이션
- **Phase 1/2/3 임상시험 설계**: 단계별 PK 파라미터 예측
- **PK 지표 산출**: Cmax, AUC, t½, Vss 등

### 3. 📋 IND Application Generator
- **FDA Form 1571 스타일 보고서 생성**
- 프로젝트 데이터 자동 매핑 (PK, 독성, 시뮬레이션 결과)
- PDF 출력 지원 (Print to PDF)
- 폼 데이터 실시간 미리보기 (iframe 기반)

### 4. 📊 프로젝트 관리
- 프로젝트 생성/관리
- 코호트별 시뮬레이션 결과 저장
- Phase별 결과 비교 대시보드

---

## 📂 Project Structure

```
PKSmart/
├── app/                        # Web Application (FastAPI)
│   ├── main.py                 # 앱 진입점
│   ├── database.py             # DB 연결 설정
│   ├── models.py               # DB 모델 정의
│   ├── routers/                # API 라우터
│   │   ├── auth.py             # 인증
│   │   ├── projects.py         # 프로젝트 관리
│   │   ├── cohorts.py          # 코호트 관리
│   │   ├── ind_agent.py        # IND Generator
│   │   └── dashboard.py        # 대시보드
│   ├── services/               # 비즈니스 로직
│   │   └── ind_generator.py    # IND 문서 생성 서비스
│   ├── static/                 # CSS, JS
│   └── templates/              # HTML (Jinja2)
│       ├── ind_generator.html  # IND Generator 페이지
│       ├── ind_fda1571.html    # FDA 1571 양식 템플릿
│       ├── project_results.html# 프로젝트 결과
│       └── cohort_detail.html  # 코호트 상세
├── pksmart/                    # AI Core (mPBPK, QSAR)
│   ├── mpbpk_engine.py         # mPBPK 시뮬레이션 엔진
│   ├── qsar_predictor.py       # QSAR 예측기
│   └── cyp2d6/                 # CYP2D6 대사 모듈
├── models/                     # 학습된 ML 모델 (.joblib)
├── data/                       # 데이터 파일
│   └── generated_ind/          # 생성된 IND 보고서
└── requirements.txt
```

---

## 🚀 웹 대시보드 실행 방법

### 📂 실행 경로
```
ClinicalTrials-main/PKSmart/
```

### 💻 실행 명령어
```bash
# 1. PKSmart 디렉토리로 이동
cd ClinicalTrials-main/PKSmart

# 2. 가상환경 생성 및 활성화
python -m venv .venv
.\.venv\Scripts\activate  # Windows
source .venv/bin/activate  # Linux/Mac

# 3. 의존성 설치
pip install -r requirements.txt

# 4. 서버 실행
uvicorn app.main:app --reload --port 8000

# 5. 브라우저에서 접속
# http://localhost:8000
```


---

## 🔧 주요 API 엔드포인트

| Method | Endpoint | 설명 |
|--------|----------|------|
| GET | `/dashboard` | 메인 대시보드 |
| GET | `/projects/{id}` | 프로젝트 상세 |
| GET | `/ind-generator` | IND Generator 페이지 |
| GET | `/ind-generator/export-fda1571` | FDA 1571 보고서 생성 |
| POST | `/api/ind/generate` | IND 문서 생성 API |
| POST | `/cohorts/create` | 코호트 생성 |
| GET | `/cohorts/{id}` | 코호트 상세 |

---

## � 데이터 흐름

```
[SMILES 입력]
      ↓
[QSAR 독성 예측] → Tox21 12개 엔드포인트
      ↓
[PK 파라미터 예측] → Cmax, AUC, t½, Vss
      ↓
[mPBPK 시뮬레이션] → 가상 코호트 생성
      ↓
[IND Generator] → FDA 1571 보고서 생성
```

---

## 🛠️ 기술 스택

| 분류 | 기술 |
|------|------|
| **Backend** | FastAPI, SQLAlchemy |
| **Frontend** | Jinja2, Tailwind CSS |
| **ML/AI** | scikit-learn, RDKit, Mordred |
| **Database** | SQLite |
| **LLM** | Google Gemini API |

---

## 📝 최근 업데이트 (2026.02.08)

- ✅ FDA Form 1571 스타일 IND 보고서 템플릿 추가
- ✅ IND Generator에 FDA 1571 미리보기 통합 (iframe)
- ✅ 폼 데이터 실시간 반영 (Applicant, Drug, PK/Tox)
- ✅ Expected Patients 필드에 실제 시뮬레이션 수 반영
- ✅ Print / Save as PDF 기능 추가

---

## � License

This project is for educational and research purposes.
