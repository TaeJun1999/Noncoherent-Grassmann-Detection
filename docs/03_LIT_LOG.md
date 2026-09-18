# 03_LIT_LOG — 검색으로 확인한 선행연구의 기록

> **append-only.** 검색·열람으로 확인한 것만 적는다. 기억으로 쓴 arXiv 번호는 없다.
> 항목 형식: ID / 날짜 / 질문 / 결론 / 목록 / 판단 / ⚠️ 한계.
> "확인 수준" 열의 뜻: **원문** = 본문을 직접 읽음 / **코드** = 공개 구현을 직접 읽음 / **초록** = 초록·서지 페이지만 / **2차** = 다른 논문의 서술·참고문헌으로만.
> `01_PLAN.md` §3(차별화 대상)과 역할이 다르다. §3은 논문에 쓸 서사이고, 이 파일은 그 근거다. §2.3 재선택 후 §3을 고칠 때 이 파일에서 옮긴다.
>
> 2026-09-18 신설: L1.

---

## L1 — `M>1` 구조화 Grassmannian constellation의 효율 demapper 존재 여부 (2026-09-17 ~ 09-18)

**질문** R6에서 `M=1` Cube-split의 그리디 demapper가 exact ML에 1 dB 미만으로 붙었다. 1차 주장 (i)안이 살려면 `M>1`에 같은 급의 싼 demapper가 없어야 한다. 있는가.

### 결론

**있다. 단 전부 hard decision 전용이다.**

- `M=2`에는 Grass-Lattice(그리디), 일반 `M`에는 Exp-Map simplified decoder와 training/PSAM 계열이 있다. 셋 다 비용이 SVD에 지배된다.
- `M>1` 구조화 codebook에서 그리디급 비용으로 **LLR을 내는 demapper는 찾지 못했다.** CP code 쪽은 고차원 복호와 soft 복호를 미해결로 명시한다.
- 같은 `M>1` codebook 위에서 싼 demapper 대 exact ML을 비교한 공개 수치는 하나다: 직교화 PSAM, `T=4, M=N=2`, 약 1 dB.
- 일반 `M`으로 가는 기하 기반 설계(Cube-split, Grass-Lattice)는 멈춰 있다. Cube-split은 부록 A에서 미해결로 남겼고, Grass-Lattice는 `M=2`까지만 나왔다.

### A. `M>1`, 구조화, 그리디급 demapper가 있는 것

| # | 논문 | 식별자 | `M` | 싼 demapper (비용) | soft | ML 대비 격차 (공개값) | 확인 수준 |
|---|---|---|---|---|---|---|---|
| A1 | Cuevas–Beltrán–Gutiérrez–Santamaría–Tuček, "Structured Multi-Antenna Grassmannian Constellations for Noncoherent Communications," IEEE SAM 2024 | arXiv 없음. IEEE Xplore 10636457. 같은 내용이 Cuevas 박사논문(EURASIP 공개) §4.3.2 | 2 전용 | `svd(Y)` → 2×2 SVD → 스칼라 근 찾기 3회 → 정규 CDF `4(T−2)`회 → 좌표별 반올림. SVD 지배 | 없음 | 없음. 저자가 ML 비용을 이유로 비교를 생략 | 박사논문 원문 + 코드 (grassbox) |
| A2 | Kammoun–Cipriano–Belfiore, "Non-Coherent Codes over the Grassmannian," IEEE TWC 6(10):3657–3667, 2007 (Exp-Map) | arXiv 없음. DOI 10.1109/TWC.2007.06059 | 일반 | `svd(Y)` → `M×M` SVD로 CS 형태 정렬 → `Ĥ = X̂ᴴY` → `Y·Ĥ⁻¹`의 아래 블록 → QAM slicing `M(T−M)`개. SVD 지배 | 없음 | 있음 (초록이 GLRT와 simplified 두 복호기 평가를 언급). 수치 미열람 | 초록 + 코드 (grassbox, `M≤2` 특수화) |
| A3 | Calabuig, "Unitary Space-Time Modulation Based on Coherent Codes," Sensors 22(23):9049, 2022 (직교화 PSAM) | arXiv 없음. DOI 10.3390/s22239049. PMC9738684 | 일반 정식화, 평가는 2 | LS 채널 추정 + Schnorr–Euchner sphere decoder. §5.1의 2단계 수신기: `T/M`개 PSAM 각각에서 후보 하나 → 후보 중 ML | 없음 | **약 1 dB.** `T=4, M=N=2`, `L=256`과 `65,536`. 비부호화 block error rate. 싼 수신기는 직교화를 모르는 mismatched 수신기 | 원문 (Europe PMC 전문) |
| A4 | Dayal–Brehler–Varanasi, "Leveraging coherent space-time codes for noncoherent communication via training," IEEE Trans. IT 50(9):2058–2080, 2004 | arXiv 없음 | 일반 | MMSE 추정 + coherent 복호. 저비용이 설계 목표 | coherent soft demapper를 붙일 수 있으나 평가 미확인 | 미확인 | 초록 |
| A5 | Ashikhmin–Calderbank, "Grassmannian Packings From Operator Reed–Muller Codes," IEEE Trans. IT 56(11):5689–5714, 2010 | arXiv 없음 | `2^k` | 이진 RM 복호 알고리즘의 수정판 (초록). order 미확인 | 미확인 | 미확인 | 초록 |

codebook 크기 (식에서 계산):
- A1: `B = 4B'(T−2)` bit, `B'`는 실수 성분당 비트. `T=4`: 8, 16, 24 bit (`B'=1,2,3`). `T=6, B'=1`: 16 bit.
- A2: `B = M(T−M)·log₂Q` bit. `M=2, T=4`, 4-QAM: 8 bit. `T=16, M=4`, 4-QAM: 96 bit.
- A3: `(L_Q)^{2M(T−M)}`점, `L_Q = √Q`. `T=4, M=2`, 4-QAM: 8 bit, 16-QAM: 16 bit. `T=16, M=2`, 16-QAM: 112 bit. 평가는 `T∈{4,8,16}`.

### B. `M>1`이지만 그리디급이 아니거나, 복호가 미해결인 것

| # | 논문 | 식별자 | 내용 | 확인 수준 |
|---|---|---|---|---|
| B1 | Gohary–Davidson, "Noncoherent MIMO Communication: Grassmannian Constellations and Efficient Detection," IEEE Trans. IT 55(3):1176–1205, 2009 | arXiv 없음. McMaster 공개본 | 비구조 codebook + `K` 크기 lookup table의 감축 탐색. `K=1024`, 30 dB에서 평균 51.2회 likelihood 평가 (ML은 1024회). 평가 `K ≤ 1024`. hard 전용. 비용이 `K`에 비례하므로 `B ≥ 20`에서 그리디급이 아니다. baseline #4의 원형 | 원문 |
| B2 | El-Azizy–Gohary–Davidson, "A BICM-IDD scheme for non-coherent MIMO communication," IEEE TWC 8(2):541–546, 2009 | arXiv 없음 | B1 기반 list demapper + 복호기 feedback으로 list 보강. **`M>1` soft 출력의 유일한 비-exhaustive 선례.** 비용은 `K`에 비례 | 초록 |
| B3 | Asano 외, "Sparse Grassmannian Design for Noncoherent Codes via Schubert Cell Decomposition" | arXiv 2601.21009 (2026-01) | exhaustive GLRT를 `O(\|X\|TN)`으로 (`M`배 절감). `2^B`에 선형. 서론이 "기하 기반 설계는 대부분 SIMO, MIMO의 성능–검출 복잡도 균형은 여전히 major challenge"라고 서술 | 원문 |
| B4 | Soleymani–Mahdavifar, "Analog Subspace Coding…," IEEE Trans. IT 68(4), 2022 (CP codes). 고차원 확장은 같은 저자의 ISIT 2021 "New packings in Grassmannian space" | arXiv 1909.07533 | 다항식 평가 기반 구조화 subspace code | 초록 |
| B5 | Riasat–Mahdavifar, "Decoding Analog Subspace Codes: Algorithms for Character-Polynomial Codes" | arXiv 2407.03606 (2024-07) | RS 복호를 이용한 최소거리·list 복호. **대상은 1차원 CP code, hard decision뿐.** 결론부: 고차원 CP code 복호는 block code로 직접 사상되지 않아 "may require developing entirely new techniques". soft-decision도 향후 과제로 명시 | 원문 |
| B6 | Pendás-Recondo 형제, "A Structured Family of Grassmannian Constellations via Geodesic Mapping…" | arXiv 2510.15070 | 최대 `4M²`점. `B = 20–24`와 무관 | 초록 |
| B7 | Lanham 외 / 일반화판, QEC 기반 noncoherent STBC | arXiv 1812.07115, 2305.07104 | 저율. ML 복호 규칙 유도 | 초록 |
| B8 | Pauli–Lampe, "Tree-search multiple-symbol differential decoding for unitary space-time modulation," IEEE TCOM 55(8):1567–1576, 2007 | arXiv 없음 | 차분 unitary ST의 다중심볼 검출에 tree search. baseline #4 선례. soft판(Pauli–Lampe–Schober, IT 52(4), 2006)은 단일 안테나 DPSK용 | 초록 |

### C. `M=1` (비교 기준)

| # | 논문 | 식별자 | 내용 | 확인 수준 |
|---|---|---|---|---|
| C1 | Ngo–Decurninge–Guillaud–Yang, Cube-split, IEEE TWC 19(3), 2020 | arXiv 1905.08745 | 그리디 복호 + LLR. **부록 A 원문: `M>1` 확장은 "not evident and left as perspective for future work".** Voronoi cell을 좌표 조건으로 쓸 수 없다는 것이 이유 | 원문 |
| C2 | Cuevas 외, Grass-Lattice, IEEE TWC 23(3), 2024 | arXiv 2209.04172 | `O(T²N)` 복호, hard 전용. greedy 대 ML 비교는 `M=1`, `T∈{2,4}`만. 결론부: ϑ₁·ϑ₃은 일반 `M`으로 유도, **ϑ₂의 MIMO 확장은 향후 과제** | 원문 |
| C3 | Beltrán–Ferizović–López-Gómez, "Measure-preserving mappings from the unit cube to some symmetric spaces," J. Approx. Theory 308, 2025 | arXiv 2303.00405 | measure-preserving 사상은 rank-one 대칭공간(구, 사영공간)까지. `G(M,T)`, `M>1`은 범위 밖 | 초록 |
| C4 | Shigenaga 외, Z-Opt | arXiv 2605.04545 (2026-05) | `G(2,1)` 전용 | 초록 |

### D. 인접 — 모델은 다르지만 reviewer가 물을 것

| # | 논문 | 식별자 | 내용 | 확인 수준 |
|---|---|---|---|---|
| D1 | Ngo–Guillaud–Decurninge–Yang–Schniter, EP 다중사용자 검출, IEEE TWC 19(9), 2020 | arXiv 1905.11152 | `Y = Σ xₖhₖᵀ + W`. soft 출력, `K, N, T`의 다항식 비용. 열이 직교하지 않아 `XᴴX ≠ I`이므로 확정된 Bingham 형태 밖이다. 하지만 "rank-`M` noncoherent soft 검출을 다항식 비용으로"는 선점돼 있다. 차별화 문장 필요 | 초록 |
| D2 | Baba–Iimori–Pradhan–Malomsoky–Ishikawa, "Deep Neural Network Based Reduced-Complexity Detector for Grassmann Constellation," IEEE VTC-Fall 2024 | arXiv 미발견 | 학습 기반 검출기. baseline #9(discriminative)의 직접 선례. `M`, 입력 표현, 비교 대상 미확인 | 2차 |

### E. 지형 교차 확인에 쓴 자료

- 서베이: Ngo–Cuevas–de Miguel Gil–Monzon Baeza–García Armada–Santamaría, "Noncoherent MIMO Communications: Theoretical Foundation, Design Approaches, and Future Challenges," arXiv 2505.23172 (2025-05). Table IV가 구조화 설계의 전체 목록이다. 2020년 이후 `M>1` 구조화 항목은 A1, B4, B7뿐이다. **원문**
- Cuevas, "Advanced Grassmannian Constellation Designs for Noncoherent MIMO Communications," 박사논문, Universidad de Cantabria, 2024 (EURASIP 공개). **원문**
- grassbox toolbox, `github.com/diegocuevasfdez/grassbox`, BSD-3, Zenodo DOI 10.5281/zenodo.14857110. Grass-Lattice(`M∈{1,2}`), Exp-Map(`M∈{1,2}`, `M=2`는 `T∈{4,6}`), Cube-split의 인코더·복호기 MATLAB 구현. **코드**

### §2.3 재선택에 넘기는 판단

1. **(i)안의 hard decision 판은 약하다.** 공개된 유일한 격차가 문턱값 1 dB이고, 그것도 mismatched 수신기 기준이다. R6에서는 10³ MAC급 좌표강하 보정이 `M=1` 격차를 0.3–0.6 dB에서 ≤ 0.1 dB로 줄였다. 같은 보정이 `M=2`에서 통하는지는 측정된 적이 없다. 외삽하지 않는다.
2. **빈 자리는 "`M>1` soft 출력"이다.** A1–A5 중 LLR을 내는 싼 demapper는 없다. B5는 미해결이라고 스스로 적었다. 남는 선례는 `K`에 비례하는 B2와 다른 모델의 D1이다. 다만 A3·A4에 coherent soft MIMO demapper를 붙이는 조합은 자명한 경쟁자이므로 baseline으로 세워야 한다. 평가된 사례는 찾지 못했다.
3. **헤드라인 `M=4`의 구조화 경쟁자는 A2, A3/A4, A5다.** A1은 `M=2`에서 멈춰 있다 (C2, C3).
4. **새 gate는 바로 구현 가능하다.** `M=2, T=4`에서 A1, A2, A3 모두 8–16 bit codebook이 나오므로 exhaustive ML이 돈다 (R6이 `B=17–20`을 돌렸다). grassbox가 참조 구현이다. `00_STATE.md` `[VERIFY]` 5번.

### ⚠️ 한계

- A3의 "약 1 dB"는 저자의 본문 서술이다. 그림은 보지 못했고 어느 block error rate 수준인지 모른다.
- A2, A4, A5, B2, D2는 본문 미열람이다 (IEEE 유료). 우선순위: A5(`M=4` 해당) > A2(격차 수치, 일반 `M` 복호기 단계) > D2.
- A2의 비용 서술은 grassbox의 `M≤2` 구현에서 읽은 것이다. 일반 `M`에서도 같은 구조일 것으로 보이나 원문으로 확인하지 않았다 `[VERIFY]`.
- 검색은 영어 문헌, 2026-09-18 시점이다. "없다"는 "찾지 못했다"의 뜻이다.
- 2026-09-17 중간 보고에 "같은 codebook에서 ML 격차를 잰 논문은 하나도 없다"고 적었다. A3 본문을 읽고 철회했다. 이 파일이 최종본이다.
