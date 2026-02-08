# FDA 1571 스타일 IND 보고서 템플릿 구현 계획

## 목표
현재 IND Generator의 출력물을 FDA Form 1571과 유사한 공식 양식 스타일로 변환

## FDA 1571 양식 분석 (3페이지 구조)

### Page 1: 기본 정보
| 필드 | 설명 | 데이터 소스 |
|------|------|-------------|
| 1. Center | CDER/CBER 선택 | 고정값 |
| 2. Name of Sponsor | 스폰서명 | `project.title` |
| 3. Date of Submission | 제출일 | 현재 날짜 |
| 4. Sponsor Address | 주소 | 입력 필드 |
| 5. Telephone Number | 전화번호 | 입력 필드 |
| 7A. IND Number | IND 번호 | 자동 생성 |
| 7B. IND Type | Commercial/Research | 체크박스 |
| 8. Name of Drug | 약물명 | `prediction.smiles` 기반 |
| 8A. Indication for Use | 적응증 | 입력 필드 |
| 9. Phase | Phase 1/2/3 | `cohort.phase` |
| 12A. Submission Type | 제출 유형 | 체크박스 |

### Page 2: 신청 내용 체크리스트
| 항목 | 21 CFR 참조 | 현재 데이터 매핑 |
|------|------------|----------------|
| 1. Form FDA 1571 | 312.23(a)(1) | ✅ 본 양식 |
| 2. Table of Contents | 312.23(a)(2) | 자동 생성 |
| 3. Introductory Statement | 312.23(a)(3) | 약물 개요 |
| 4. General Investigative Plan | 312.23(a)(3) | 연구 계획 |
| 5. Investigator's Brochure | 312.23(a)(5) | PK 요약 |
| 6. Protocol | 312.23(a)(6) | 프로토콜 정보 |
| **7. CMC Data** | 312.23(a)(7) | 제조 정보 |
| **8. Pharmacology/Toxicology** | 312.23(a)(8) | **IND Generator PK 데이터** |
| 9. Previous Human Experience | 312.23(a)(9) | 기존 시험 정보 |
| 10. Additional Information | 312.23(a)(10) | 기타 정보 |

### Page 3: 서명 정보
| 필드 | 설명 |
|------|------|
| 19. Sponsor Representative | 스폰서 대표자 정보 |
| 20-21. Contact | 전화/팩스 |
| 22. Address | 주소 |
| 23. Email | 이메일 |
| 24. Signature Date | 서명일 |
| 28-29. Signatures | 서명란 |

---

## 구현 방안

### Option A: HTML/CSS 기반 인쇄용 템플릿 (권장)
- 새 템플릿 `ind_fda1571.html` 생성
- CSS `@media print` 스타일로 PDF 출력 최적화
- 브라우저 인쇄 기능으로 PDF 저장

### Option B: Python PDF 생성 (ReportLab/WeasyPrint)
- 백엔드에서 직접 PDF 생성
- 더 정교한 레이아웃 제어 가능
- 추가 라이브러리 필요

---

## 제안: Option A (HTML/CSS 기반)

### 수정/생성 파일
1. **[NEW] `app/templates/ind_fda1571.html`** - FDA 1571 스타일 템플릿
2. **[MODIFY] `app/routers/ind_agent.py`** - `/ind-generator/fda1571` 라우트 추가
3. **[NEW] `app/static/css/fda1571.css`** - 인쇄용 스타일

### UI 통합
- IND Generator 페이지에 "Export as FDA 1571" 버튼 추가
- 클릭 시 새 탭에서 FDA 1571 스타일 페이지 열림
- Ctrl+P로 PDF 저장 가능

---

## 다음 단계
1. 사용자 승인 후 구현 시작
2. HTML 템플릿 생성
3. CSS 스타일링
4. 라우터 연결
5. 테스트 및 검증
