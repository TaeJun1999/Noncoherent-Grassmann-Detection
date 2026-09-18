# 00_STATE — 현재 상태

> **충돌하면 이 파일이 이긴다.** 마일스톤마다 전체 교체한다.
> 마지막 갱신: 2026-09-18 (R6 (e) 측정, [LIT] L1 직후)

---

## 이번 교체에서 바뀐 것

- "지금의 gate" 근거 둘째·셋째 불릿 교체. 둘째는 R6 (e) 반영, 셋째는 `03_LIT_LOG.md` L1 반영. **"`M>1`에는 싼 demapper가 문헌에 없다"는 문장은 틀렸다. 삭제했다.**
- "확정" 표의 R6 행에 R6 (e) 결과 추가.
- 잔여 `[VERIFY]` 3번을 닫음. 5번(`M=2` gate), 6번(미열람 본문) 추가.
- "1주차·R6가 바꾼 설계" 표의 1차 주장 행을 [LIT]·R6 (e) 기준으로 다시 씀.
- "다음 작업" 갱신. 2번·4번 완료 처리, `M=2` gate [EXP] 추가.
- "실행 환경" 절 신설 (저장소, 서버 실행 절차).
- 복잡도 체크리스트 마지막 항목 교체.
- "폐기" 목록은 손대지 않았다. 저자 확인 없이 폐기 항목을 늘리지 않는다.
- `01_PLAN.md`에 **stale 표시만** 넣었다 (본문 0줄 변경). §2.3, §2.4, §3, §4에 `⚠️ STALE` 상자. 내용 교체는 §2.3 재선택 이후.
- "작업 분담" 절 신설. Claude / GPT / Git의 역할과 파일 배치.

---

## 지금의 gate

**R6 판정: 구조화 codebook의 효율 demapper와 exact ML의 격차는 1 dB 미만이다. gate 실패.**

`01_PLAN.md` §2.3 무대 선택으로 복귀한다. §2.3 재선택이 끝나기 전에는 네트워크 구현·학습 작업을 시작하지 않는다.

근거 (`02_RESULTS_LOG.md` R6, Cube-split `M=1`, exhaustive ML이 ground truth):

- hard decision(SER): 원논문 그리디 복호기가 0.3–0.6 dB, 비용 10³ MAC급 좌표강하 보정을 붙이면 ≤ 0.1 dB. `T ∈ {4,8}`, `B ∈ {8,14,17,20}`, `N ∈ {1,2,4}`, seed 두 개에서 같다.
- soft 출력(LDPC BLER): `T=4`에서 ≤ 0.3 dB. `T=8`에서 원논문 식 (26) 리스트 demapper 대비 +0.7 dB이며, 리스트 크기 η를 16→32로 늘려도 0.05 dB만 준다. R6 (e): 리스트 중심을 near-ML 점(D3)으로 옮겨도 +0.46 dB가 남는다. 중심 문제가 아니다. 격차는 후보 집합 밖의 posterior mass와 함께 움직인다 (잡힌 mass 0.786 / 0.813 / 0.892에서 BLER 0.2533 / 0.1500 / 0.0667, exact 0.0367). `T`개 cell 리스트의 합집합은 exact 비용의 약 1/46로 +0.21 dB까지 접근한다. dB 값은 한 점 BLER의 환산값이다.
- `M>1`에도 구조화 codebook의 싼 hard demapper가 있다 (`03_LIT_LOG.md` L1). `M=2` 전용 Grass-Lattice(SAM 2024), 일반 `M`의 Exp-Map simplified decoder(TWC 2007), 직교화 PSAM + LS/sphere decoder(Sensors 2022). 셋 다 SVD 지배 비용이고 **soft 출력이 없다.** 공개된 ML 대비 격차는 하나뿐이다: 직교화 PSAM, `T=4, M=N=2`, 약 1 dB (비부호화, 저자 서술값, mismatched 수신기). Cube-split 원논문 부록 A가 `M>1` 확장을 미해결로 남긴 것은 사실이나, 그것이 "문헌에 없다"를 뜻하지는 않는다.

R5의 3.09 dB는 대용 shortcut이 만든 숫자였다. 더 이상 인용하지 않는다.

---

## 확정 — 다시 유도하지 않는다

| 항목 | 내용 | 근거 |
|---|---|---|
| Posterior | `p(X\|G) ∝ p(X)·exp(κ tr(XᴴGX))`, `κ = ρ/(σ²(ρ+σ²))` | R1 |
| 형태 | Haar prior이면 complex matrix Bingham(`κG`) | R1 |
| 정규화 상수 | `₁F̃₁(M;T;κG)`, `G`의 고유값에만 의존 | R2 |
| 충분통계량 | `G`는 zero-mean Gaussian family에서 충분. 공간 상관에서도 유지 | R1 |
| 충분성 붕괴 | Rician · 위상잡음 · 송신 비선형 | R1 |
| Eigenframe | `A = VᴴX`, 조건 입력은 `λ` 하나. U(T)-equivariance는 구조적 | R1 |
| MACG base | `c = 2κ/T` (측정 `c* ≈ 1.8κ/T`) | R3 |
| Gold standard | MACG 제안 independence MH, acceptance ~0.5 | R2 |
| **구조화 codebook의 싼 demapper** | `M=1` Cube-split에서 hard decision은 near-ML (≤ 0.6 dB, 보정 시 ≤ 0.1 dB). soft 격차는 `T=8`에서 ~0.7 dB, η로도 중심 이동으로도 닫히지 않음. `T`-cell 합집합(후보 4352개)은 +0.2 dB까지 접근 | R6, R6 (e) |

## 폐기 — 되살리지 않는다

- "관측 = forward process의 한 시점이므로 guidance 없이 reverse 시작"
- `N = M` clean special case
- `σ_t`로 인덱싱된 물리 경로를 conditional path로 사용
- singular Wishart 밀도 — 이 문제에서는 필요하지 않다
- `span(Y)`만으로 조건화
- **"구조화 codebook에서 싼 demapper가 정확도를 잃는다"는 전제** (R5는 대용품, R6가 부정)
- **D3s** — 다른 좌표를 hard decision에 고정한 soft 확장. 과신으로 LDPC를 망친다 (R6)

## 잔여 `[VERIFY]`

1. `01_PLAN.md` §1.5-3의 Duistermaat–Heckman 측도 형태 (불변성 자체는 검증됨)
2. `01_PLAN.md` §2.5 두 번째 화살표: LLR 섭동 → BLER regret
3. ~~`T=8` soft 격차 0.7 dB의 성격~~ → **닫힘 (R6 e).** 중심 문제가 아니라 놓친 posterior mass 문제다
4. **`T ≥ 16` codebook posterior의 gold standard** — exhaustive 불가(`B ≥ 34`). 필요한지 여부는 §2.3 결정에 달림
5. **`M=2, T=4` gate** — Grass-Lattice, Exp-Map, 직교화 PSAM 세 codebook에서 싼 demapper 대 exact ML. SER과 LDPC BLER. 코드 없음. (i)안 채택의 전제
6. **[LIT] 미열람 본문** — Ashikhmin–Calderbank 2010(`M=2^k`, 헤드라인 `M=4` 해당), Kammoun 2007(GLRT 대 simplified 수치), Baba 2024(baseline #9 선례). `03_LIT_LOG.md` L1의 ⚠️ 참조

---

## 1주차·R6가 바꾼 설계

| 항목 | 이전 | 현재 | 근거 |
|---|---|---|---|
| 헤드라인 설정 | `T=8, M=2, B=16` | `T=16–32, M=4, B=20–24` (§2.3 재선택까지 잠정) | R4 |
| 네트워크 폭 | 미정 | `d ≤ 64` | R4 |
| NFE | "few-step 바람직" | ≤ 4가 목표 | R4 |
| base 분포 | Haar | MACG(`c=2κ/T`). 고 SNR ESS 5–20× | R3 |
| 1차 주장 문구 | "동일 비용에서 정확도 격차를 회수" | **근거 상실. §2.3 재선택 대기.** 남는 안: (i) `M>1` 한정 — hard decision 경쟁자가 셋 있고 공개 격차가 문턱값(약 1 dB)이다. soft 출력은 문헌에 비어 있다; (ii) mismatch(2차)를 1차로 승격; (iii) `T≥16` soft 출력 무대 — R6 (e)로 폐기되지 않았으나 `T=8, M=1`에서는 남은 자리가 +0.2 dB뿐이다. gold standard 없이는 주장 불가 | R6, R6 (e), L1 |

---

## 다음 작업 (chat 단위)

1. ~~[STATE] R6 반영~~ — 완료 (2026-09-17).
2. ~~[LIT] `M>1` 구조화 demapper~~ — 완료 (`03_LIT_LOG.md` L1). 답: 있다, 단 hard 전용. 잔여는 `[VERIFY]` 6번.
3. **[DERIV] §2.3 재선택** — 입력이 다 모였다: R6, R6 (e), L1, R4 break-even. (i)/(ii)/(iii) 중 결정. 결정 chat이지 실험 chat이 아니다. 산출물은 `01_PLAN.md` §2.3–2.4 교체 문안. 세 안 각각에 "죽는 조건"을 적는다.
4. ~~[EXP] `p26_center_test.py 1.5 300 16`~~ — 완료 (R6 e). (iii)안 폐기 조건 미충족.
5. **[EXP] `M=2, T=4` gate** — `[VERIFY]` 5번. 3번 앞에 둘지 뒤에 둘지는 저자 결정. 앞에 두면 (i)안 판정에 숫자가 생긴다.
6. **[DERIV]** codebook posterior gold standard (`T≥16`) — 3에서 (iii)을 남긴 경우에만.

3이 끝나기 전에는 figure·네트워크·학습 작업을 시작하지 않는다. 5번은 측정이므로 이 제한에 걸리지 않는다.

---

## 작업 분담 (2026-09-18 확립)

세 곳을 쓴다. **Git이 원본이고, 나머지 둘은 마일스톤 사본이다.** 사본의 머리말에 그 시점을 적는다. 사본과 Git이 다르면 Git의 최신 커밋이 새 것이다.

| | 역할 | 가진 것 |
|---|---|---|
| **Claude project** | 판단·방향·gate. 문헌 확인. 문서 교체본 생성. 코드 작성과 짧은 검산 | `00_STATE.md`, `01_PLAN.md`, `02_RESULTS_LOG.md`, `03_LIT_LOG.md` |
| **GPT project** | 긴 유도, 증명, 대안 계산, 코드 초안 | 위 문서 4개 + `code/cubesplit_gap.py` + 그때 다루는 스크립트 |
| **GitHub** | 원본. 코드, 원시 로그, 문서, GPT 산출물 | 전부 |

규약:

- **GPT 산출물은 chat에 붙여넣지 않는다.** `docs/derivations/D<n>_*.md`로 커밋하고 commit hash만 Claude에 준다. Claude가 clone해서 읽는다.
- GPT가 낸 것은 제안이다. 확정이 아니다. 수식은 Claude나 서버에서 검산한 뒤에야 `00_STATE.md` "확정"에 올라간다.
- GPT는 숫자와 arXiv 번호를 새로 만들지 않는다. `02_RESULTS_LOG.md`와 `03_LIT_LOG.md`에 있는 것만 인용한다.
- 문헌 확인은 Claude가 한다. 한 곳에 몰아야 "검색하지 않은 논문을 인용하지 않는다"가 지켜진다.

---

## 실행 환경 (2026-09-18 확립)

- 저장소: `github.com/TaeJun1999/Noncoherent-Grassmann-Detection` (**public 유지**, `main`). 코드는 `code/`, 결과는 `results/<ID>_<설정>/`, 문서는 `docs/`.
- 긴 실행은 저자 서버에서 Claude Code가 돌린다. 절차: pull → `sanity.py` → 본 실행(배경 + 폴링) → 원시 로그만 커밋 → `REPORT` 블록(`code_commit`, `result_commit` 포함).
- Claude는 `result_commit` hash로 고정해 clone한 뒤 `stdout.log` 원문에서만 숫자를 옮긴다. 코드 동일성은 sha256으로 확인한다.
- 서버의 PATH `python3`에는 numpy가 없다. `~/miniforge3/envs/torch/bin/python`을 절대경로로 쓴다 (저장소 `CLAUDE.md` 규칙).
- `lists_*.npy`는 `.gitignore`에 있다. 커밋하지 않는다.

---

## figure 먼저 (§2.3 결정 이후)

초안에 가장 먼저 필요한 것은 그림이다. 두 개를 먼저 그린다. 단, 무대가 바뀌면 2번 그림의 서사도 바뀐다.

1. **문제 설정 그림** — block-fading, `G`로의 축약, Grassmann 위 posterior
2. **핵심 아이디어 그림** — MACG base에서 출발한 flow가 mismatch에서만 크게 움직이는 모습. "이동량 = 모델 오차" 서사를 한 장에

---

## 복잡도 주장 체크리스트

논문에 복잡도 수치를 쓸 때마다 확인한다.

- [ ] `n_sc`와 심볼 길이를 명시했는가 (없으면 검증 불가능한 수치다)
- [ ] 우리 학습 비용과 baseline의 setup 비용에 같은 회계를 적용했는가
- [ ] 무엇을 세고 무엇을 제외했는지, 그 이유를 적었는가
- [ ] NFE와 latency와 복소 곱셈 수를 모두 보고했는가
- [ ] 경쟁 demapper의 리스트 크기 `η`와 그 비용을 함께 적었는가 (R6: `2·B·η` 후보 × `T·N` MAC)
- [ ] `M>1` 경쟁 demapper(Grass-Lattice, Exp-Map, 직교화 PSAM)의 비용(SVD 지배)과 soft 출력 부재를 출처(`03_LIT_LOG.md` L1)와 함께 명시했는가
