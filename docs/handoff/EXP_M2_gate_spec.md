# EXP 사양 — `M=2, T=4` gate (R7)

> `00_STATE.md` `[VERIFY]` 5번의 실행 사양. 작성 2026-09-18.
> 숫자의 출처는 `docs/02_RESULTS_LOG.md`, 문헌 근거는 `docs/03_LIT_LOG.md` L1.
> 이 문서는 "무엇을 어떤 정의로 재는가"만 정한다. 결론과 해석은 여기 쓰지 않는다.

---

## 1. 질문

`M>1` 구조화 Grassmannian codebook에서, 문헌의 싼 demapper와 exact ML의 격차는 얼마인가.

R6이 `M=1` Cube-split에서 답한 것과 같은 질문이다. 비교 가능해야 하므로 지표와 절차를 R6에
맞춘다.

**완료 판정.** 두 codebook 각각에서 (a) 비부호 SER 격차 dB, (b) 부호화 BLER 격차 dB.
각 격차는 수신기별로 분리해 보고한다.

**이 숫자가 정하는 것.** `01_PLAN.md` §2.3 (i)안의 생사. 격차가 R6과 같은 규모(1 dB 미만)면
(i)안은 `M=1`과 같은 결말이다. 격차가 크게 남으면 (i)안에 근거가 생긴다.

---

## 2. 범위

**이번에 재는 것: A1 Grass-Lattice(`M=2`), A2 Exp-Map(`M=2`).** 둘 다 `03_LIT_LOG.md` L1의
grassbox에 참조 구현이 있고, `T=4`에서 8–16 bit codebook이 나와 exhaustive ML이 돈다.

**이번에 재지 않는 것: A3 직교화 PSAM.** 공개 구현이 없어 논문에서 새로 구현해야 한다.
게다가 A3의 공개 격차 약 1 dB는 직교화를 모르는 mismatched 수신기 기준이다. R6에서 싼
수신기를 조금만 고쳐도 격차가 0.3–0.6 dB에서 ≤ 0.1 dB로 줄었다. 그 교훈을 A1·A2에 먼저
적용하는 편이 정보가 많다. A3는 A1·A2 결과를 보고 결정한다.

---

## 3. 신호 모델과 정규화 — 가장 먼저 고정할 것

R6의 `M=1` 규약은 `Y = √(ρT)·x hᵀ + Z`, `‖x‖=1`, `h ~ CN(0,I_N)`, `Z ~ CN(0,1)`이다.
`M`으로 확장할 때 **송신 전력을 채널 사용당 `ρ`로 유지한다.**

```
Y = √(ρT/M) · X H + Z,    X ∈ C^{T×M}, XᴴX = I_M,  H ~ CN(0, I) (M×N),  Z ~ CN(0,1)
Σ = (ρT/M)·XXᴴ + I_T
κ = (ρT/M) / (1 + ρT/M)
ML 지표 = tr(Xᴴ G X),   G = YYᴴ
```

`M=1`에서 R6 식으로 환원되는지 코드로 확인한다 (`sanity` 항목 1).

**SNR 축의 뜻을 표에 반드시 적는다.** `ρ`는 채널 사용당·수신 안테나당 SNR이다. 블록 SNR은
`ρT`다. `M=1` 결과와 dB를 나란히 놓을 때 이 정규화 때문에 `10log₁₀M` 만큼의 상수 차이가
생긴다. 표 제목에 규약을 적지 않으면 나중에 읽을 수 없다.

---

## 4. Codebook

| | A1 Grass-Lattice `M=2` | A2 Exp-Map `M=2` |
|---|---|---|
| 출처 | Cuevas 외, IEEE SAM 2024. 같은 내용이 Cuevas 박사논문 §4.3.2 | Kammoun–Cipriano–Belfiore, IEEE TWC 6(10), 2007 |
| 참조 구현 | grassbox `GrassLatticeEncoding.m` / `GrassLatticeDecoding.m` | grassbox `ExpMapEncoding.m` / `ExpMapDecoding.m` |
| 크기 | `B = 4B'(T−2)`. `T=4`: `B'=1` → 8 bit, `B'=2` → 16 bit | `B = M(T−M)·log₂Q`. `T=4, M=2`: 4-QAM → 8 bit, 16-QAM → 16 bit |
| 싼 복호기 | SVD → 2×2 SVD → 스칼라 근 3회 → 정규 CDF `4(T−2)`회 → 좌표별 반올림 | SVD → `M×M` SVD로 CS 정렬 → `Ĥ = X̂ᴴY` → `Y·Ĥ⁻¹`의 아래 블록 → QAM slicing `M(T−M)`개 |

둘 다 비용이 SVD에 지배된다. 둘 다 soft 출력이 없다 (L1).

`T=4`에서 `B ∈ {8, 16}` 두 크기를 모두 돌린다. `B=16`이면 `K=65536`으로 R6의 `K=131072`보다
작다. exhaustive ML은 문제없다.

비트 라벨은 참조 구현의 Gray 규약을 그대로 쓴다. 라벨이 다르면 BLER이 달라진다.

---

## 5. 수신기 — R6과 1:1로 대응시킨다

| 이름 | 정의 | R6 대응 |
|---|---|---|
| ML | 전수 `argmax_k tr(X_kᴴ G X_k)` | ML |
| G1 | 참조 구현의 그리디 복호기 그대로 | D1 |
| G2 | G1 출력에서 시작해 격자·QAM 좌표를 하나씩 exact 지표로 좌표강하 (최대 3 sweep) | D3 |
| exact LLR | 전 codebook posterior 비트 marginal, ±30 clip | exact LLR |
| L1(η) | 비트별 chordal 최근접 `η`개 리스트(사전계산), **G1 결정 중심**, 리스트 위 log-sum | P26(η) |
| L2(η) | 같은 리스트, **G2 결정 중심** | R6 (e)의 P26-D3 |

`M>1`에는 문헌의 싼 soft demapper가 없다 (L1). 따라서 L1(η)와 L2(η)는 우리가 세우는
경쟁자다. 원논문이 준 것이 아니다. **논문에 쓸 때 이 사실을 반드시 명시한다.** 구성은
Cube-split 식 (26)과 같은 방식이며, R6 (e)에서 중심 이동이 격차를 절반 가까이 줄인 것이
확인됐으므로 L2를 빼지 않는다.

`η ∈ {8, 16}`. 리스트 크기를 키운 효과는 R6에서 이미 작았다 (0.05 dB). 크기보다 중심이
중요하다는 것이 R6 (e)의 결론이다.

**cell 합집합의 대응물은 이번에 만들지 않는다.** Cube-split의 `T`개 cell은 그 codebook의
구조다. A1·A2에 같은 구조가 있는지는 확인되지 않았다. 있으면 다음 회차에 넣는다.

---

## 6. 재사용과 신규 작성

`code/cubesplit_gap.py`에서 **그대로 쓰는 것** (codebook·`M`과 무관):

```
exact_llr, lse, snr_at, gram, LDPC, make_ldpc, gf2_rref, channel의 잡음 생성 부분
```

**`M=1` 가정이 박혀 있어 일반화가 필요한 것** (네 군데뿐이다):

| 함수 | 현재 | 필요한 형태 |
|---|---|---|
| `channel` | `X (n,T)`, `H (n,N)`, `√(ρT)` | `X (n,T,M)`, `H (n,M,N)`, `√(ρT/M)` |
| `kappa_of` | `ρT/(1+ρT)` | `(ρT/M)/(1+ρT/M)` |
| `codebook_features` | `C (T,K)`에서 `c_kᴴ G c_k` | `X_k (T,M)`에서 `tr(X_kᴴ G X_k) = ⟨G, X_kX_kᴴ⟩`. 같은 Hermitian 내적 트릭이 그대로 성립한다. 특징 차원은 여전히 `T²` |
| `precompute_lists_sorted` | 근접도 `\|c_iᴴ c_j\|²` | `‖X_iᴴ X_j‖_F²` (chordal) |

`ml_metrics_fast`, `exact_llr_chunked`, `gvec`는 `codebook_features`가 위 형태로 바뀌면
수정 없이 동작한다. `Phi`의 열 정의만 맞으면 된다.

**새로 쓰는 것**: A1·A2의 인코더, 라벨, G1 복호기, G2 좌표강하.

파일은 `code/m2_gate.py` (라이브러리), `code/sanity_m2.py` (검증), `code/run_m2_ser.py`,
`code/run_m2_bler.py`로 나눈다. R6의 파일 분할과 같은 구조다.

---

## 7. 구현 검증 (`sanity_m2.py`) — 본 실행 전에 반드시 통과

| # | 항목 | 통과 기준 |
|---|---|---|
| 1 | `M=1`로 환원. 새 `channel`·`kappa_of`·`codebook_features`에 `M=1` Cube-split codebook을 넣고 R6 함수와 대조 | 지표 상대오차 < 1e-6, 같은 seed에서 SER 동일 |
| 2 | codebook이 Grassmann 위에 있는가 | 모든 `k`에서 `‖X_kᴴX_k − I₂‖_F < 1e-10` |
| 3 | 비트 라벨 왕복 | 1000개 무작위 인덱스에서 일치 |
| 4 | 무잡음 복원율 (G1, G2) | 2000 블록에서 1.0000 |
| 5 | G2 지표 ≤ ML 최대 지표 | True |
| 6 | 참조 구현 대조 (§8) | 같은 입력에서 같은 codeword 인덱스 |
| 7 | exact ML f32 대 f64 argmax 일치 | ≥ 0.999 |
| 8 | 리스트 사전계산의 `η`번째 경계 동률 비율 | 보고만 한다. 안정 정렬로 결정론화 |
| 9 | LDPC 구성 | R6과 같은 rate 0.502 확인 |

4번이 1.0000이 아니면 인코더나 복호기가 틀린 것이다. 그 상태로 SER을 재지 않는다.

---

## 8. 참조 구현 대조

grassbox는 MATLAB이다. Octave로 돌려 참조 벡터를 만든다.

- grassbox를 **commit hash로 고정**해 받고, 그 hash를 결과에 적는다.
- 고정 seed로 인덱스 200개에 대해 codeword 행렬을, 그리고 고정 `Y` 200개에 대해 복호 결과
  인덱스를 `.npz`로 저장한다.
- Python 포트가 같은 입력에서 같은 값을 내는지 확인한다. 허용오차: codeword `1e-10`,
  복호 인덱스는 완전 일치.

**Octave를 쓸 수 없으면** 이 항목을 건너뛰고 §7의 1–5, 7–9로 대신한다. 그 경우
`⚠️ 한계`에 "참조 구현과 대조하지 못했다"를 반드시 적는다. 대조 없이 "원논문 복호기와
같다"고 쓰지 않는다.

---

## 9. 스윕 설정

R6과 같은 규약이다.

**비부호 SER**: 점당 4000 블록, seed 11. 재현으로 seed 23에서 8000 블록 한 번.
보고는 SER 1e-1과 1e-2에서의 SNR(dB), log-linear 보간. 1e-3은 오차가 커서 인용하지 않는다.

**부호화 BLER**: LDPC (3,6)-regular Gallager, rate 약 0.502, sum-product 50회, 무작위 비트
인터리버, seed 5. `n`은 `B`와 6의 공배수 (`B=8` → `n=1008`, `B=16` → `n=1008`). 점당 최소
300–400 부호어, 블록 오류 30개 이상까지 최대 1500.

**비용**: 블록당 복소 MAC을 해석적 계수로 적고, 1코어 wall time도 함께 잰다. SVD 비용을
명시한다. R6 (d)와 같은 계수를 쓴다.

---

## 10. 하지 말 것

- 참조 구현의 수식을 기억으로 재구성하기. `.m` 파일과 논문을 읽는다.
- `M=1` 결과에서 `M=2` 격차를 외삽하기.
- 검증 4번이 1.0000이 아닌 상태로 스윕 돌리기.
- 우리가 만든 L1(η)·L2(η)를 "원논문의 demapper"라고 부르기.
- A3(PSAM)를 이번 회차에 끼워 넣기.
- 결과를 보고 사양을 고치기. 고쳐야 하면 고친 사실과 이유를 결과에 적는다.
