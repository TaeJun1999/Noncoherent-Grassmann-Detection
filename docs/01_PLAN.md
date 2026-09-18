# 01_PLAN — 물리 정합 Grassmann posterior inference (noncoherent MIMO)

> 연구의 **안정 본체**: 정식화, 주장 구조, 선행연구, baseline, 설계 결정.
> 유도가 바뀌거나 설계 결정이 뒤집힐 때만 교체한다.
> 현재 상태는 `00_STATE.md`, 숫자는 `02_RESULTS_LOG.md`가 유일한 출처다.
>
> 계보: handoff v1 (2026-09-16) → v2 (Bingham 유도) → v3 (1주차 수치 검증) → 01_PLAN
>
> ---
>
> **⚠️ STALE 표시 (2026-09-18 추가). 본문은 고치지 않았다.**
>
> R6 (gate 실패), R6 (e), `03_LIT_LOG.md` L1이 나온 뒤 이 파일의 일부가 낡았다.
> 낡은 절 앞에 `⚠️ STALE` 상자를 달았다. **내용은 손대지 않았다.** §2.3 재선택([DERIV])이
> 끝나면 §2.3, §2.4, §3, §4를 한 번에 교체한다. 그 전까지 이 파일에서 다음을 인용하지 않는다:
>
> | 절 | 낡은 이유 | 대신 볼 곳 |
> |---|---|---|
> | §2.3 무대 선택 | "구조화 codebook을 주 무대로 채택"의 전제가 R6에서 무너짐 | `00_STATE.md` "지금의 gate" |
> | §2.4 1차 주장 | "동일 비용에서 정확도 격차를 회수" 근거 상실. `§8.1–8.2` 참조는 존재하지 않는 절을 가리킴 | `00_STATE.md` "1주차·R6가 바꾼 설계" |
> | §3 선행연구 표 | `M>1` 구조화 demapper 항목이 빠져 있음 | `03_LIT_LOG.md` L1 |
> | §4 baseline #4, #5 | R6에서 측정 완료. L1이 경쟁자 셋을 추가 | `02_RESULTS_LOG.md` R6, `03_LIT_LOG.md` L1 |
>
> §1(정식화), §5(설계), §6(검증 위생), §7(사실 확인)은 유효하다. §1은 R1–R3로 검증됐다.

---

## 0. 목적

1순위 주제 하나만 파고든다. 목표는 system model, 학습 objective, 실험 설계를 논문 수준으로 확정하는 것.

---

## 1. 정식화 (확정)

### 1.1 시스템 모델

Rayleigh block-fading noncoherent MIMO.

```
Y = X H + W
```

| 기호 | 정의 |
|---|---|
| `X ∈ G(M,T)` | T×M orthonormal, 송신 subspace |
| `H` | M×N 채널, 열이 i.i.d. `CN(0, ρI_M)` |
| `Y` | T×N 수신 신호 |
| `W` | T×N AWGN, 원소 `CN(0, σ²)` |

`Y`의 열은 i.i.d. `CN(0, Σ)`, `Σ = ρXXᴴ + σ²I_T`.

### 1.2 충분통계량

```
p(Y|X) = (π^T det Σ)^(−N) · exp( −tr(Σ⁻¹ YYᴴ) )
```

`Y`가 `G = YYᴴ`를 통해서만 들어간다. 따라서 **`G`는 `X`에 대한 충분통계량**이다.

이 성질은 `N`과 `T`의 대소와 무관하고, `H`에 공간 상관이 있어도 (`y_n ~ CN(0, XRXᴴ+σ²I)`) 유지된다. **깨지는 지점은 정확히 non-zero mean(Rician), 위상잡음, 송신 비선형이다.** §5.6의 조건화 ablation이 여기서 나온다.

`span(Y)`만 남기면 특이값 가중을 버린다. 저 SNR에서는 subspace가 아니라 수신 신호의 스펙트럼이 지배적이라는 분석이 이미 나와 있다 (estimator–correlator 동작과 일치; arXiv 2505.08432). §3 표에 반영.

### 1.3 Posterior closed form ← 이 주제의 새 출발점

`det Σ = (ρ+σ²)^M σ^{2(T−M)}`는 `X`와 무관하다. Woodbury로

```
Σ⁻¹ = σ⁻²( I_T − β XXᴴ ),   β = ρ/(ρ+σ²)
tr(Σ⁻¹G) = σ⁻²·tr G − (β/σ²)·tr(XᴴGX)
```

따라서

```
p(X|G) ∝ p(X) · exp( κ · tr(XᴴGX) ),   κ = ρ / ( σ²(ρ+σ²) )
```

**Haar prior에서 이것은 `G(M,T)` 위의 complex matrix Bingham, parameter `κG`이다.** GLRT metric `tr(XᴴYYᴴX)`는 이 밀도의 지수부 그 자체다 — 즉 GLRT는 이 posterior의 MAP이다.

### 1.4 정규화 상수

`∫ exp(κ tr(XᴴGX)) [dX] = ₁F̃₁(M; T; κG)` — 복소 행렬변수 hypergeometric 함수이며 `G`의 고유값에만 의존한다. 이 함수는 Jacobi ensemble 밀도에 대한 linear spectral statistic `tr(YG)`의 생성함수로 적분 표현을 가진다 (arXiv 1508.01358). 수치 계산은 Koev–Edelman 알고리즘. saddlepoint 근사는 실제 차원에서 정밀도가 부족하다는 보고가 있으니 gold standard 용도로는 쓰지 말 것.

> §9(v1의 §8) 4순위 "Jacobi ensemble 좌표"가 이 적분 표현과 같은 대상이다. 두 주제가 예상보다 가깝다.

### 1.5 Eigenframe 축약과 불변성

`G = VΛVᴴ`, `A = VᴴX`로 두면

```
tr(XᴴGX) = tr(AᴴΛA) = Σ_t λ_t · ‖row_t(A)‖²
```

세 가지가 따라온다.

1. **U(T)-equivariance가 공짜.** `V`로 회전 → `A` 생성 → `X = VA`. 학습으로 equivariance를 얻을 필요가 없다.
2. **조건 입력이 `λ ∈ R^T₊` 하나로 붕괴한다** (+ `σ², ρ, N`). T×T 행렬을 네트워크에 넣을 이유가 없다.
3. 밀도가 `A ↦ Dᴴ A Q` (`D` 대각 유니터리, `Q ∈ U(M)`)에 불변이다. 즉 `A`에 대한 의존이 moment map `r(A) = (‖row_t(A)‖²)_t ∈ Δ(M,T)` (hypersimplex)를 통해서만 일어난다. **matched 문제는 실질적으로 T차원 문제다.**

3번은 논문에서 방어가 아니라 공격 재료로 쓴다: "언제 생성 추론이 불필요한가"를 정확히 특정하는 절이 된다. non-Haar prior에서 이 축약이 깨지며, **깨짐의 정도가 곧 novelty의 크기다.**

`[VERIFY]` 3번의 Duistermaat–Heckman 측도 형태(hypersimplex 위 조각별 다항식)는 미확인. 불변성 자체는 자명하므로 축약 논증에는 영향 없음.

### 1.6 폐기된 정식화 (기록용)

- **"관측 = forward process의 한 시점, guidance 없이 reverse 시작"** — v1에서 이미 폐기. `span(Y)`가 충분통계량이 아니라서.
- **`N = M` clean special case** — v2에서 추가 폐기. `N=M`이어도 `span(Y)`는 여전히 고유값을 버린다. 충분통계량은 `(span(Y), 양의 고유값 M개)`다. 게다가 forward process 해석이 성립하려면 `{MACG(ρXXᴴ+σ²I)}_σ`가 semigroup이어야 하는데 그럴 근거가 없다. **이론 섹션에서도 쓰지 말 것.**
- **`σ_t`로 인덱싱된 물리 경로를 conditional path로 사용** — 같은 이유로 폐기. marginal이 맞지 않는다.

---

## 2. 범위와 주장 (재구성)

### 2.1 범위

고정 rank · single-user · noncoherent soft decoding. 다중 사용자 · rank 변화 · semantic JSCC는 전부 후속 논문.

### 2.2 matched-Haar 정리가 강제하는 것

§1.3의 귀결을 정면으로 받아야 한다.

**Haar prior + matched Rayleigh에서 baseline 3(parametric Bingham)은 근사가 아니라 exact다.** 더구나 Bingham은 unimodal이다 — critical point는 `G`의 고유벡터 M개를 고르는 `C(T,M)`개이고, eigengap이 있는 한 global max는 top-M 하나뿐이다. 다중 모드가 없다.

즉 이 regime에서 **"diffusion이 parametric보다 정확하다"는 주장은 성립할 수 없다.** v1 §4가 이를 "가장 위험한 반론"으로 적었는데, 반론 수준이 아니라 증명 가능한 사실이다.

### 2.3 무대 선택

> **⚠️ STALE.** 아래 표의 "주 무대로 채택"은 R6 이전 판단이다. R6에서 구조화 codebook의 싼
> demapper가 exact ML에 1 dB 미만으로 붙어 gate가 실패했다. 세 무대의 재선택이 진행 중이다
> ([DERIV], `00_STATE.md` "다음 작업" 3번). 아래 판정 열을 인용하지 않는다.

논문이 설 수 있는 곳은 셋이다.

| 무대 | posterior가 Bingham이 아닌 이유 | 판정 |
|---|---|---|
| **구조화된 discrete codebook prior** | `p(X)`가 `2^B`개 점에 지지 → 연속 Bingham이 아니라 조합적 marginalization 문제 | **주 무대로 채택.** §2.4 1차 주장과 정합 |
| **Mismatch** (Rician, phase noise, PA 비선형, site-specific) | likelihood가 깨지고 `G`의 충분성도 깨짐 | **2차 주장으로 채택.** site-specific을 headline으로 |
| `N < M` | `rank(G)=N<M` → Bingham parameter가 rank-deficient → maximizer가 양차원 집합 | 진짜 비단봉이지만 capacity-suboptimal이라 설득력이 약함. **eigengap→0이라는 연속 knob으로 흡수** |

### 2.4 주장

> **⚠️ STALE.** 1차 주장의 근거가 R6에서 무너졌다. 아래 "측정된 위치 (v3, §8.1–8.2)" 불릿의
> 참조는 이 파일에 존재하지 않는 절을 가리킨다 (v3 시절의 잔재). 헤드라인 `T=16–32, M=4,
> B=20–24`는 R4 비용 분석에서 나온 값이라 유효하지만, 그 위에 얹힌 주장 문구는 §2.3 재선택
> 대기 중이다. 2차 주장(mismatch)과 "Matched case의 새 역할"은 영향을 받지 않았다.

**1차 주장.** 물리 정합 생성 추론이 noncoherent soft decoding의 정확도–계산량 관계를 개선한다.

- matched model의 exact MAP보다 낮은 Bayes error를 목표로 하지 **않는다**. 더 적은 계산으로 그 성능에 접근하는 것이 목표.
- 무대: 구조화된 대형 constellation, coded soft output 등 exact 탐색이 비싼 영역.
- **측정된 위치 (v3, §8.1–8.2).** 제안 수신기는 두 극단 사이를 메운다: 구조 없는 exhaustive보다 `B ≥ 15`에서 싸고, 구조화 효율 demapper보다는 약 10⁴배 비싸다. 따라서 주장은 "가장 싸다"가 아니라 **"동일 정확도를 훨씬 싸게, 또는 동일 비용에서 정확도 격차를 회수한다"**로 쓴다. 헤드라인 설정은 `T=16–32, M=4, B=20–24`.
- **잔여 전제 하나가 아직 미측정이다. §8.2의 ⚠️를 먼저 볼 것.**

**2차 주장.** likelihood가 깨지는 mismatch 환경에서 학습된 모델이 강건하다.

- "likelihood를 아는데 왜 diffusion인가"에 대한 답은 이제 §2.2가 정확히 규정한다: **likelihood를 정확히 알면 Bingham으로 끝난다. 그래서 무대는 likelihood를 모르는 곳이다.**
- 이때 "정확한 커널"은 closed form이 아니라 **학습에 쓴 simulator(QuaDRiGa / Sionna RT) 기준**임을 반드시 명시.

**Matched case의 새 역할.** 성능 무대가 아니라 **검증 무대**다. unnormalized posterior가 정확하고 정규화 상수도 계산 가능하므로, posterior 근사 오차를 exact하게 측정할 수 있다. 대부분의 생성모델 논문이 못 하는 검증이다. §6의 자산으로 쓴다.

### 2.5 이론 (논문 내부, 별도 프로젝트 아님)

- I-MMSE 일반화를 노리지 말 것. (v1 유지)
- **목표 교체.** "일반 posterior 근사 오차 → decoding regret"은 몇 주 안에 닫히지 않는다. 대신:

```
missing posterior mass  →  per-bit LLR 오차  →  BLER regret
```

리스트가 놓친 집합 `S`에 대해 LLR 오차가 `log(1 + p(S)/p(list))` 꼴로 제어되고, `p(S)`는 SNIS의 χ²-divergence와 리스트 크기 `L`로 bound된다. 작은 `B`에서는 brute force로 직접 측정까지 가능하므로 **이론과 실험이 같은 축에 놓인다.**

`[VERIFY]` 두 번째 화살표(LLR 섭동 → BLER regret)의 형태. mismatched decoding / GMI 쪽 표준 논증으로 닫힐 것으로 보이나 미확인. 몇 주 안에 안 닫히면 부록으로 축소.

---

## 3. 반드시 차별화해야 할 선행연구

> **⚠️ STALE (불완전).** 아래 10개는 유효하다. 다만 `M>1` 구조화 constellation의 효율
> demapper가 빠져 있다. Grass-Lattice(SAM 2024), Exp-Map(TWC 2007), 직교화 PSAM(Sensors 2022),
> El-Azizy BICM-IDD(TWC 2009), Ashikhmin–Calderbank(IT 2010)를 넣어야 한다. 확인 수준과
> 서지 정보는 `03_LIT_LOG.md` L1에 있다. 옮겨 적지 말고 그 파일을 인용한다.

v1의 6개 + v2 추가 4개.

| 논문 | 내용 | 차별화 / 용도 |
|---|---|---|
| **Land-then-transport** (arXiv 2601.07512) | 물리 채널을 continuous-time flow에 통합. AWGN Gaussian path + closed-form teacher velocity. 추론 시 유효 잡음 분산으로 ODE 시작 시각 설정 | "채널을 generative path에 넣었다"는 Euclidean에서 끝난 이야기. 우리 차별점: 비동기라 등화 불가 / unknown이 manifold 위 / 커널이 non-Gaussian / **posterior가 Bingham으로 닫히고 그 축약이 언제 깨지는지를 규정** |
| **GRAM-DIFF** (2602.15187) | Gram matrix를 reverse process의 guidance 항으로 주입 | 우리는 Gram을 guidance가 아니라 exact conditioning으로. 게다가 **eigenframe에서 `λ`로 붕괴**시킨다 |
| **DiffATS** (2605.09275) | Procrustes 정렬로 gauge 제거 후 Euclidean diffusion | 정렬 후 Euclidean은 필수 baseline (#5). gauge 제거만으로는 novelty 없음 |
| RSGM / Riemannian Flow Matching / Grassmann Manifold Flows | manifold 위 생성모델 일반론 | 기성품. heat kernel 회피도 기여가 아님. **RFM 채택 근거로만 인용** |
| CSI feedback 계열 (2605.08628, 2606.06846) | 생성 prior 기반 CSI feedback | feedback은 붐빔 → 이 논문은 detection |
| Besson, *MMSD of a subspace* (arXiv 1101.3462) | Bingham prior + MCMC exact Bayesian subspace estimation | **exact baseline (#2)** |
| **Low-complexity noncoherent massive MIMO detection** (arXiv 2505.08432) | 저 SNR에서 subspace가 아니라 스펙트럼이 지배적임을 KL 분석으로 제시 | **§1.2 관측의 외부 근거.** v1에 없던 인용 |
| **Hoff (2009), matrix Bingham–vMF simulation** | Stiefel 위 rejection/Gibbs sampler | **gold-standard sampler로 채택** |
| **Kent–Ganeiber–Mardia (2013)** | ACG envelope 기반 Bingham rejection sampling; MACG 밀도와 샘플링 구성 | §5.2 base 분포의 근거 |
| **Cube-split** (Ngo–Decurninge–Guillaud–Yang, TWC 2020) / **Gohary–Davidson** (IT 2009) | 구조화 Grassmannian constellation과 효율적 검출 | **baseline #8의 실체.** 이게 없으면 계산량 주장이 검증되지 않음 |

---

## 4. Baseline 세트 (10개, 하나도 빼지 말 것)

> **⚠️ STALE (부분).** #4와 #5는 "없으면 1차 주장의 전제가 미검증"이라고 적혀 있다. R6에서
> `M=1` Cube-split에 대해 측정했고 결과는 gate 실패다 (`02_RESULTS_LOG.md` R6, R6 (e)).
> #5의 `M>1` 대응물은 L1이 셋을 찾았다. 나머지 8개는 유효하다.

**고전 / exact**

1. GLRT / exact ML — 작은 constellation에서 exact posterior와 정확도 대조
2. Besson류 Bingham prior + MCMC (MMSD)
3. **Parametric MACG / Bingham posterior** — matched case에서는 exact임을 명시하고, 그 사실을 논문이 먼저 말할 것
4. **트리 / sphere-search 기반 noncoherent soft demapper** `[v2 추가]`
5. **codebook 구조를 쓰는 효율 demapper** (Cube-split 계열 per-bit demapper) `[v2 추가]` — 4와 5가 없으면 1차 주장의 전제가 미검증

**학습 기반**

6. Euclidean diffusion
7. 정렬 후 Euclidean diffusion (DiffATS 방식)
8. Deterministic Grassmann network
9. **Discriminative amortized LLR network** `[v2 추가]` — `(λ, V)`에서 LLR을 직접 회귀. **parametric 다음으로 위험한 baseline.** 이게 비기면 생성 기계 전체가 정당화되지 않는다

**제안 방법과 그 ablation**

10. 제안 intrinsic 모델 (§5)
    - 10a. 학습 제거 ablation: `MACG(I+cG)` proposal + exact rescoring만

**비교 규칙**

- 반드시 동일 계산량(NFE / latency) 축에서 비교
- oracle evaluation(정답 채널 보고 후보 선택) 절대 금지
- 모든 방법에 **같은 회계**를 적용한다. 우리 학습 비용을 빼면서 baseline의 setup 비용을 넣지 말 것. 무엇을 세고 무엇을 제외했는지, 그 이유와 함께 명시
- baseline을 그룹으로 소개하고(고전 / 학습), 각 그룹에 대해 **이득의 출처를 따로 귀속**할 것 (고전 대비: 학습 능력 / 학습 baseline 대비: 기하 구조와 exact rescoring)

---

## 5. 확정된 설계

### 5.1 조건화 네트워크

§1.5의 eigenframe 축약을 그대로 구현한다.

```
G = VΛVᴴ  →  A = VᴴX 공간에서 생성  →  X = VA
조건 입력: λ ∈ R^T₊  (+ σ², ρ, N)
```

속도장 parametrization:

```
v_θ = (I − AAᴴ) · F_θ(P, λ) · A,    P = AAᴴ
```

- `U(M)` gauge equivariance와 horizontality가 구조적으로 보장된다 (`P`는 gauge 불변, 좌측 사영이 horizontality를 강제).
- `F_θ`는 T개 토큰 위의 transformer. 토큰 feature `(λ_t, P_tt)`, pairwise feature `|P_ts|` + 위상.
- `D(T)`-equivariance 조건 `F(DᴴPD, λ) = Dᴴ F(P,λ) D`는 위상을 `P_ts/|P_ts|`로 실어 나르면 만족한다.
- `T ≤ 32`에서 `T²` attention은 비용이 아니다.

**`N < T`는 확률론적 문제가 아니다.** §1.2의 충분성 논증이 Gaussian 밀도에서 직접 나오므로 singular Wishart 밀도를 쓸 일이 없다. `p(G)`를 모델링하지 않기 때문이다. 남는 것은 구현 이슈뿐 — 인코더가 0 고유값 `T−N`개를 가진 rank-deficient spectrum을 처리해야 한다. **v1 §5-2와 §7-2는 삭제한다.**

### 5.2 학습 objective — RFM, base는 Haar가 아니라 MACG

**Riemannian conditional flow matching 채택.** score matching은 배제한다 — RSGM 계열은 manifold 위 heat kernel을 Sturm–Liouville 전개나 Varadhan 점근으로 근사해야 하는 반면, RFM은 conditional flow를 벡터장의 적분곡선으로 정의하고 manifold 거리로 매칭한다. Grassmann은 exp/log가 SVD로 닫혀 있다.

- conditional path: geodesic interpolant. 목표 속도 `u_t = log_{A_t}(A_1) / (1−t)`.
- **base 분포를 uniform이 아니라 `MACG(I + cG)`로 잡는다.** 샘플링은 `Z ~ CN(0,I)` (T×M), `X = orth((I+cG)^{1/2} Z)`.

base 변경이 가져오는 것:

1. parametric posterior가 baseline이 아니라 **모델의 구성요소**가 된다. "unimodal이면 parametric으로 충분"이라는 반론이 구조적으로 무력화된다.
2. matched case에서 flow가 near-identity가 되므로 **few-step이 설계상 자연스럽다.** 4–8 NFE가 현실적 목표.
3. mismatch에서의 이동량이 곧 "물리 모델이 틀린 정도"의 측정치다. 그대로 논문 figure가 된다.

**`c`는 측정으로 확정됐다** (`code/scale_and_cost.py`, part A).

```
c* ≈ a · κ/T,   a = 1.8 (중앙값), 범위 1.06–6.0
기본값 c = 2κ/T,  그 다음 (T,M,SNR)당 ESS line search 한 번
```

작은 `κ`에서는 `a ≈ 1.1`로 깨끗하고, `κ`가 커질수록 `a`가 흩어진다. 유도 근거: MACG 지수부를 1차로 전개하면 `≈ T·c·tr(XᴴGX)`이므로 `Tc ≈ κ`.

**부수 결과 — base 선택의 독립적 근거.** Haar base 대비 튜닝된 MACG base의 ESS 이득:

| κ | Haar ESS/n | MACG(c*) ESS/n | 이득 |
|---|---|---|---|
| 1.0 | 0.86–0.96 | 0.999 | 미미 |
| 4.0 | 0.04–0.46 | 0.76–0.98 | **5–20×** |

즉 **고 SNR에서 base 변경만으로 sample efficiency가 한 자릿수 개선된다.** flow의 기여와 무관하게 성립하므로, ablation #10a(학습 없는 proposal)가 이미 의미 있는 baseline임을 뜻한다.

### 5.3 Soft output — sampler를 posterior가 아니라 proposal로 쓴다

```
샘플 {X⁽ˡ⁾} ~ q_θ(·|G)
  → codebook nearest (chordal)
  → dedupe
  → exact metric κ·tr(Xᵢᴴ G Xᵢ)로 재채점
  → max-log LLR
```

재채점이 들어가면 `q_θ`는 정확할 필요가 없고 **coverage만 좋으면 된다.** 정확도는 exact metric이 보장하고, 계산량은 리스트 크기가 결정한다. 1차 주장과 정확히 맞물리는 구조다.

- 비트별 빈 리스트를 막기 위해 hard-decision codeword를 강제 포함시킨다.
- 누락 가설에 대한 LLR 처리(floor / clipping) 규칙을 명시하고 baseline에도 동일 적용.
- 외부 코드는 5G NR LDPC 또는 polar. **BLER을 보고하되, 채널 코딩 없는 수치를 "end-to-end 성능"이라 부르지 말 것.**

### 5.4 Few-step 추론

base가 정보를 담고 있으므로 distillation 이전에 이미 저 NFE에서 동작해야 한다. 순서:

1. 먼저 NFE 스윕으로 Pareto 곡선을 그린다 (정확도 × NFE).
2. 부족하면 Riemannian consistency / progressive distillation. exp/log가 닫혀 있어 그대로 이식된다.
3. 보고는 NFE와 **고정 GPU에서의 latency, 그리고 복소 곱셈 수** 둘 다. 계산량 축이 1차 주장이므로 하나만으로는 부족하다.

### 5.5 조건화 ablation — `G` vs full `Y`

`G`는 zero-mean Gaussian family 전체에 대해 충분하지만 Rician · 위상잡음 · 비선형에서는 lossy하다. 따라서:

| 조건 입력 | 의미 |
|---|---|
| `(U_Y, s)` = 좌특이벡터 + 특이값 (`G`와 동치) | 물리 충분통계량. 불변성 최대 |
| full `Y` | 충분성이 깨진 regime에서만 추가 정보 |

이 ablation이 **"언제 물리 충분통계량이 도움이고 언제 족쇄인가"**에 답한다. 리뷰어가 물어보기 전에 우리가 먼저 답하는 항목.

### 5.6 Mismatch 실험 설계

| impairment | 유효성 | 함정 / 조건 |
|---|---|---|
| spatial correlation | **약함** | `y_n ~ CN(0, XRXᴴ+σ²I)`도 Gaussian이라 `R`을 알면 여전히 closed form. **반드시 `R` unknown / 소표본 추정 조건으로** |
| Rician, unknown LoS 위상 | 강함 | non-central → `G`의 충분성 붕괴 |
| per-symbol phase noise | 강함 | block-constant 가정 파괴, likelihood가 적분 형태 |
| PA 비선형 (Rapp / Saleh) | 강함 | `Y ≠ XH` |
| **site-specific (Sionna RT / DeepMIMO)** | **headline** | 유효 prior가 non-Haar → §1.5 축약이 깨짐 → 진짜 다봉 |

각 impairment마다 **matched baseline을 정직하게 준다**: nominal `σ², ρ`를 쓴 closed-form Bingham/ACG, 그리고 가능하면 robustified 고전 방법(유효 잡음 분산 추정 등). baseline의 한계는 결과 해석 전에 본문에서 미리 밝혀 둔다.

---

## 6. 검증 위생

- **matched case gold standard를 먼저 만든다.** Hoff Gibbs sampler + Koev–Edelman 정규화 상수. 제안 모델의 posterior 근사 오차를 χ² / TV / ESS로 **정확히** 측정한다. 이 자산이 §2.4 "검증 무대" 주장의 실체다.
- 10단 baseline으로 geometry 효과 / diffusion 효과 / exact rescoring 효과를 분리한다. 특히 **#9(discriminative)와 #10a(학습 없는 proposal)가 분리의 핵심축**이다.
- diffusion이 필요한 조건(non-Haar prior, 작은 eigengap, mismatch)을 먼저 특정하고, **parametric으로 충분한 경계를 논문이 먼저 명시한다.** §2.2가 그 경계를 이미 정확히 규정하므로, 방어가 아니라 기여로 쓸 수 있다.
- 평가 지표는 chordal distance / NMSE가 아니라 **BLER, soft-decoding 성능, NFE / latency / 복소 곱셈 수**.
- 합성 모델(rank · eigengap · 모드 수 통제)로 핵심 가설 검증 → site-specific(DeepMIMO / Sionna RT) 확장.
- 플롯이 무언가를 평균했다면 무엇에 대한 평균인지 항상 명시한다 (전 SNR? 전 시나리오?).

---

## 7. 사실 확인 상태

**검증됨 (웹 검색)**

- §3 선행연구 전부 (v1 6개 + v2 4개)
- MACG / GLRT 형태
- 복소 matrix Bingham 정규화 상수가 `₁F̃₁`이고 Jacobi ensemble 적분 표현을 가진다는 것
- Hoff sampler, Kent ACG envelope, Koev–Edelman의 존재

**수치 검증 완료 (v3, `code/verify_posterior.py`)**

| 항목 | 결과 |
|---|---|
| `det Σ = (ρ+σ²)^M σ^{2(T−M)}` | 상대오차 1.4e-15 |
| Woodbury `Σ⁻¹` | 절대오차 1.3e-15 |
| `log p(Y|X) − κ·tr(XᴴGX)`가 `X`에 무관 | 500회 산포 4.3e-14 (스케일 81) |
| **codebook exact posterior = Bingham weights** | **50블록 최대 TV = 4.5e-15** |
| eigenframe 축약 `tr(XᴴGX) = ⟨λ, rownorms⟩` | 8.9e-15 |
| `A ↦ DᴴAQ` 불변성 | 5.3e-15 |
| 공간 상관에서 `G` 충분성 유지 | `loglik(Y) − loglik(YU) = 0` (정확히) |
| Rician에서 `G` 충분성 붕괴 | 차이 43 nats |

**§1.3은 닫혔다.** 문서 전체가 이 유도에 얹혀 있었고, 이제 근거가 있다.

**해결됨 (v2 유도, v3에서 불필요 확인)**

- ~~Wishart 조건화의 정확한 정규화 상수~~ → §1.4. `M=1` closed form과 IS 추정이 일치 (상대오차 1e-4–2e-2)
- ~~singular Wishart 케이스~~ → §5.1. 문제 자체가 성립하지 않음
- ~~`N = M` 특수 케이스의 엄밀성~~ → §1.6. 폐기
- ~~복소 Grassmann에서 최적 `c`~~ → §5.2. `c ≈ 2κ/T`

**미검증 `[VERIFY]` (잔여 2건)**

- §1.5-3의 Duistermaat–Heckman 측도 형태 (불변성 자체는 검증됨)
- §2.5 두 번째 화살표: LLR 섭동 → BLER regret
  - 2026-09-18 보강: R6 (e)가 한 SNR 점에서 네 점을 준다. 후보 집합이 잡은 posterior mass
    0.786 / 0.813 / 0.892 / 1에 대해 BLER 0.2533 / 0.1500 / 0.0667 / 0.0367. 축이 맞다는
    증거이지 법칙이 아니다. `00_STATE.md`의 잔여 `[VERIFY]` 목록이 최신이다.

---

## 8. 실행 계획

실행 계획과 측정 결과는 이 파일에 두지 않는다. 자주 바뀌기 때문이다.

- 현재 critical path, gate, 주차 계획 → **`00_STATE.md`**
- 비용 프로파일링·검증·실험 수치 전체 → **`02_RESULTS_LOG.md`**
- 재현 코드 → **`code/`**

수치가 필요하면 기억이 아니라 `02_RESULTS_LOG.md`에서 가져온다.


---

## 9. 보류한 후속 (참고용, 이번 논문에서는 다루지 않음)

| 순위 | 주제 | 한 줄 메모 |
|---|---|---|
| 2 | Rate-weighted flag diffusion | rate-regret 가중치 `α(λᵢ−λⱼ)/(1+αλᵢ)`가 `λᵢ→λⱼ`에서 퇴화 = rank 전이 신호. must-beat baseline은 "covariance 직접 생성 후 rank 결정" |
| 3 | Rate-weighted decision-aware probing | posterior tangent 불확실성 × rate 가중치. MU interference 방향으로 좁힘. 2×2 ablation 필수 |
| 4 | Rare-outage + Jacobi ensemble 좌표 | principal angle이 failure 좌표. **§1.4의 `₁F̃₁` 적분 표현과 같은 대상** — 1순위와 이론적으로 연결됨 |
| 5 | Noncoherent URA | receiver-only, set posterior |
| 6 | ISAC 결정 층 | Wishart 조건화와 정합. SBL / gridless를 이겨야 함 |
| 7 | UL–DL bridge | 식별성 질문이 중심. 실측 cross-band 데이터 있으면 상승 |
| 8 | Pilot-free semantic JSCC | 짧은 coherence block regime 정량화가 선행 조건 |
| 9 | Hardware-realizable 생성 | 후속 응용 |

---

## 10. Venue 계획

- **첫 논문**: 위 범위 그대로 → TCOM / TWC / JSAC
- **두 번째**: 두 번째 task(ISAC few-snapshot 등)를 더한 method 논문 → NeurIPS / ICLR

용어 규칙 (원고 전체 고정)

- "utilize" 금지 → use / exploit
- "superior" 금지 → 정량화된 이득으로 대체
- "optimal" / "near-optimal"은 증명된 경우에만
- 채널 코딩 없는 수치는 "symbol detection performance"이지 "end-to-end performance"가 아니다
- 다루지 않은 능력(예: fixed-point 구현)은 언급하지 않는다
