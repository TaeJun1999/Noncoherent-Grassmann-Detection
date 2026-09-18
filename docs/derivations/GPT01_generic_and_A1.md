# [DERIV] GPT01 — 공통 층 M 일반화와 A1 Grass-Lattice 이식

**상태: 제안 / 서버 검산 전에는 확정 아님.** 작성일 2026-09-18.
기준 저장소는 `TaeJun1999/Noncoherent-Grassmann-Detection`, 읽은 기준 커밋은
`5fd114a21750354b370525e51d8fe5170e933c35`이다. 산출물은 이 문서 하나이며,
`00_STATE.md`, 결과 로그, 실행 라이브러리를 수정하지 않는다.

## 1. 질문과 범위

[EXP 사양](../handoff/EXP_M2_gate_spec.md) §3·4·6에 따라 공통 함수 네 개를
`M>=1`로 일반화하고, A1의 `M=2` 인코더·Gray 라벨·그리디 복호기를 배치 구현한다.
A2, G2, 새로운 soft demapper, SER/BLER 스윕, 성능 추정은 포함하지 않는다.
이 문서의 숫자 리터럴은 차원·입력 fixture·수치 허용오차·사양의 검산 기준이며,
실험 성능이나 예상 격차가 아니다.

## 2. 전제와 출처

[00_STATE.md의 확정 표](../00_STATE.md), R1의 posterior와 충분통계량은 인용만 한다.
아래 κ 계산은 그 식에 이번 전력 정규화를 대입하는 계산이며 posterior 재유도가 아니다.
문헌 근거는 [03_LIT_LOG.md](../03_LIT_LOG.md) L1의 **A1 및 E의 grassbox 항목**뿐이다.
논문 본문을 새로 열람했다고 주장하지 않는다. 아래 사상의 직접 근거는 첨부 실행 소스다.
새 문헌이 필요한 주장에는 `[LIT] 확인 필요`를 붙여야 하며 여기서는 새 논문을 추가하지 않는다.

코드 줄 주석의 약어는 다음과 같다. 줄 번호는 첨부 파일의 빈 줄·주석도 포함하는 1-기반 번호다.
`E=GrassLatticeEncoding.m`, `D=GrassLatticeDecoding.m`, `H=thnt.m`,
`BG=Bin2Gray.m`, `GB=Gray2Bin.m`, `F=FindAlphaOpt.m`.
`NEW`는 대응 실행문이 없는 어댑터·검증·일반화 또는 명시한 수치적 재표현이다.
`NEW(D124)`는 D124의 방정식을 다른 수치 방식으로 푼 신규 줄이지 직역이라는 뜻이 아니다.
모든 Python 실행문에 출처 또는 `NEW`가 있다. 라이선스 고지와 빈 줄은 실행문이 아니다.

**격자 값은 사양에 아직 없다.** F91, F99–100은
`lattice[j]=alpha+j(1-2alpha)/(Q-1)`를 사용하지만 gate의 alpha를 고정하지 않는다.
따라서 생성자에는 `lattice`가 필수이며 midpoint 격자나 최적 alpha를 임의로 넣지 않았다.
뒤의 alpha fixture는 F58의 탐색 범위에 있는 입력을 검산에만 사용한다.
실제 gate 전에는 저자가 격자 또는 alpha를 설정 파일에 고정하고 같은 검산을 다시 실행해야 한다.
E/D의 입력 주석은 lattice 길이를 좌표 수처럼 적지만, 실제 E55–57과 F91,100에서
lattice 길이는 **Q**다. 포트는 이 실행부를 따른다.

### R6와 맞추는 인터페이스

| 대상 | 계약 |
|---|---|
| `channel(rng,X,N,rho)` | 기존 `(n,T)`와 신규 `(n,T,M)`를 받아 `(n,T,N)` 반환. 인자 순서 그대로 |
| `kappa_of(rho,T,M=1)` | 기존 선행 인자·스칼라 의미 유지. `rho`는 선형 SNR |
| `codebook_features(C)` | 기존 `C=(T,K)` 또는 신규 `C=(T,M,K)` → `Phi=(K,T²)`, float32 |
| `precompute_lists_sorted(cs,C,bits,eta_max,chunk=256,topm=2048,log=None)` | 기존 순서 그대로. `(lists,nfall)` 반환, lists는 `(K,B,2,eta_max)` int32 |
| `GrassLattice(T,B0,lattice)` | B0는 사양의 B′. `from_index`, `index`, `bits`, `from_bits`, `symbols`, `codebook` 명칭을 R6에 맞춤 |
| A1의 `codebook()` | `(C,bits)` 반환. C의 codeword 축은 R6처럼 마지막, bits는 `(K,B)` uint8 |
| `decode_G1(gl,G)` / `decode_G1_Y(gl,Y)` | 각각 `(n,T,T)` / `(n,T,N)`에서 **0-기반 `(n,)` int64 인덱스** 반환 |

A1에는 Cube-split cell이 없으므로 라벨 메서드에 가짜 cell을 만들지 않는다.
`from_index(k)`는 `(n,nd)` digits, `index(digits)`는 `(n,)`,
`symbols(digits)`는 `(n,T,2)`다. 이는 신규 A1 표현의 차이이며 기존 R6 메서드는 바꾸지 않는다.
`N=T`이면 배열 모양만으로 Y와 G를 구분할 수 없어 입력 종류를 명시하는 두 API로 나눴다.
`exact_llr`, `lse`, `snr_at`, `gram`, `LDPC`, `make_ldpc`, `gf2_rref`,
`gvec`, `ml_metrics_fast`, `exact_llr_chunked`는 **원본 import**만 한다.
기존 `ml_metrics`, `metric_of`, `soft_P26_multi`의 rank-one C 접근까지 일반화됐다고 주장하지 않는다.
이번 코드에 새로운 soft 복호기를 넣지 않았다.

## 3. 공통 층 유도

### 3.1 전력 정규화와 κ — STATE/R1에 대입

`a=rho*T/M`라 두면 사양의 모델은
\[
Y=\sqrt a\,XH+Z,\quad X^HX=I_M,\quad H_{mn},Z_{tn}\sim\mathcal{CN}(0,1).
\]
\[
\mathbb E\|\sqrt aXH\|_F^2=a\,\mathbb E\|H\|_F^2=aMN=\rho TN.
\]
따라서 채널 사용당·수신 안테나당 신호 대 잡음 에너지 비는 rho이며,
수신 안테나당 블록 신호 에너지는 rho*T다.
[STATE/R1]의 확정 식 `κ=γ/[σ²(γ+σ²)]`에서 이 모델의 **γ=a**, **σ²=1**을 대입하면
\[
\boxed{\kappa=\frac{\rho T/M}{1+\rho T/M}}.
\]
`M=1`이면 `a=rho*T`, `κ=rho*T/(1+rho*T)`이고
`XH=x hᵀ`이므로 R6 신호식으로 환원된다. 코드도 난수 호출의 순서와 M=1 곱셈 순서를 보존한다.
SNR을 dB로 기록할 때는 `10log10(rho)`라고 명시해야 한다.
`10log10(a)=10log10(rho*T)-10log10(M)`는 축/계수 변환일 뿐,
M 변화에 따른 SER 곡선의 일정한 이동이나 성능 격차를 예측하는 식이 아니다.

### 3.2 trace 항등식과 gvec의 부호

`P_k=X_kX_k^H`라 두고 Hermitian 내적을
`<A,B>=Re tr(A^H B)`로 정의하면, Hermitian G에 대해
\[
\operatorname{tr}(X_k^HGX_k)=\operatorname{tr}(GP_k)=\langle G,P_k\rangle.
\]
첫 등호는 trace 순환성이고 둘째는 `G=G^H`와 trace의 실수성이다.
이 항등식 자체는 X가 semiunitary일 필요도 없다.
다만 R6 `codebook_features`의 내부 P는 이 P_k가 아니라
\[
\widetilde P_{ij,k}=\sum_m\overline{X_{im,k}}X_{jm,k}
=\overline{(P_k)_{ij}}=(P_k)_{ji}
\]
이다. 그러므로 실제로 저장할 특징은 R6와 동일하게
\[
\Phi_k=(\widetilde P_{ii,k};\ \Re\widetilde P_{ij,k};\ \Im\widetilde P_{ij,k})_{i<j},
\quad
g(G)=(G_{ii};\ 2\Re G_{ij};\ -2\Im G_{ij})_{i<j}.
\]
\[
\Phi_kg(G)^T=\sum_iG_{ii}P_{ii,k}
 +2\Re\sum_{i<j}G_{ij}\widetilde P_{ij,k}
 =\operatorname{tr}(GP_k).
\]
특징수는 `T+2*T(T-1)/2=T²`다. **P_k의 상삼각 허수부를 그대로 넣으면 부호가 틀린다.**
`np.triu_indices(T,1)`의 순서와 diagonal → real upper → imaginary upper의 순서를 유지한다.
따라서 R6 `gvec`·빠른 ML·chunked LLR은 수정하지 않는다.

### 3.3 chordal 리스트와 안정 정렬

모든 codeword의 rank가 같은 M이고 semiunitary이면
\[
s_{ij}=\|X_i^HX_j\|_F^2=\operatorname{tr}(P_iP_j),\qquad
 d_c^2=M-s_{ij}=\tfrac12\|P_i-P_j\|_F^2.
\]
따라서 s를 내림차순으로 정렬한다. `M=1`에서는 `s=|c_i^Hc_j|²`다.
정규화 거리 `d_c²/M`를 쓰더라도 이웃 순서는 같지만 구현에 별도 재스케일링은 넣지 않는다.
자기 자신도 해당 bit-value 리스트에 포함한다.

R6의 **상위 topm 풀 → bit/value별 stable sort → 부족한 행의 전체 후보 fallback**을 유지한다.
다만 기존 `argpartition` 출력 순서는 동률에서 인덱스 순서가 아니다.
뒤쪽 stable sort만으로는 topm 경계에서 빠진 동점 후보를 복구할 수도 없다.
이에 다음 두 줄기의 보완을 `NEW`로 표시했다: 풀의 인덱스를 먼저 오름차순으로 정렬하고,
풀 경계를 가로지르는 동률이 eta번째 선택에 걸릴 때만 전체 행으로 fallback한다.
계산된 동일 점수에서는 작은 codeword 인덱스가 먼저 온다.
이 동률 규칙 때문에 기존 R6 리스트와 동점 선택이 달라질 수 있으며 이를 M=1 수식 불일치로
취급하지 않는다. 비동점의 거리 순서는 같다. float32 계산으로 갈라진 수학적 동률까지
동점으로 합치지 않는다. tolerance 기반 임의 그룹화는 하지 않는다.

`nfall`은 R6의 실제 증가 방식대로 **(중심,bit,value) fallback 행 수**다.
고유 중심의 개수와 다르며, 이번에는 경계 동률 fallback도 포함한다.
fallback의 `for rr in bad_rows`는 Boolean 행 선택과 batched argsort로 바꾸었다.
루프는 codebook chunk·bit·bit-value뿐이며 수신 블록별 Python 루프는 없다.

## 4. A1의 각 MATLAB 단계를 수식으로 읽기

이 절의 모든 좌표 설명은 M=2, `d=T-2`, `nd=4d`, `Q=2^B0`, `B=nd*B0`다.
A1 포트는 `T>=4`를 받는다. `T=3`이면 원본의 `thnt(T-3,...)`가 양의 차수가 아니므로 지원하지 않는다.

### 4.1 E53–57, E87–102: 좌표 선택과 Gaussian 역변환

MATLAB symbol `s_j`와 Python digit `d_j`의 관계는 `s_j=d_j+1`이다.
`x_tilde_j=lattice[d_j]`이고, 좌표 순서는
\[
(\Re r_1,\Im r_1,\ldots,\Re r_d,\Im r_d,
  \Re s_1,\Im s_1,\ldots,\Re s_d,\Im s_d)
\]
에 대응하는 **확률 좌표** 순서다. 확률 좌표 x에 대해 Gaussian 성분은
`Phi_N^{-1}(x)/sqrt(2)`다. 여기서 Phi_N은 표준 실수 정규 CDF다.
E96–99의 표준편차는 1/sqrt(2)이므로 복소 r,s의 각 성분은 실수부·허수부 분산이 각각 1/2이다.
`ndtri`를 쓰며 추가 sqrt(2)를 곱하지 않는다.

### 4.2 H34–35: thnt의 방정식, 유일성, 역함수

원본의 `theta=thnt(n,t)`는 [0,1]에서 다음 식을 푼다.
\[
F_n(\theta)=(n+1)\theta^{2n}-n\theta^{2n+2}
=1-e^{-t^2}\sum_{j=0}^{n-1}\frac{t^{2j}}{j!}=:P(n,t^2).
\]
오른쪽은 정규화 lower incomplete gamma다. `z=theta²`라 놓으면
\[
F_n(\theta)=(n+1)z^n-nz^{n+1}=I_z(n,2),
\quad \frac{d}{dz}I_z(n,2)=n(n+1)z^{n-1}(1-z)>0.
\]
따라서 내부 구간에서 근은 유일하고
\[
\theta=\sqrt{I^{-1}_{P(n,t^2)}(n,2)},\qquad
 t=\sqrt{P^{-1}(n,I_{\theta^2}(n,2))}.
\]
이는 H34의 유한합에서 직접 얻은 재표현이지 새 문헌의 식이 아니다.
원본 fzero의 반복 경로를 복제하는 대신, 기존 의존성 SciPy의 벡터화된 beta/gamma
특수함수로 **동일 방정식의 근**을 구한다. 큰 확률에서는 complement를 써서 `1-p` 상쇄를 피한다.
NumPy가 배치·행렬 계산을 담당하며 별도 패키지나 scalar Python root loop를 추가하지 않는다.
MATLAB과 마지막 비트까지 같다는 주장은 아니며 참조 대조는 별도다.

영점에서 `theta/t → [(n+1)!]^(-1/(2n))`이고 역비는 그 역수다.
이는 `P(n,t²)~t^(2n)/n!`와 `I_(theta²)(n,2)~(n+1)theta^(2n)`에서 따른다.
`_radial`은 이 극한으로 0/0을 없애며, `np.where` 양쪽에서 먼저 나눗셈을 평가하지 않는다.

### 4.3 E107–118: 직교 분해와 두 radial 사상

\[
u=\frac{r^Hs}{\|r\|^2}r,\qquad \bar s=s-u,\qquad r^H\bar s=0,
\]
\[
p=\frac{r}{\|r\|}\operatorname{thnt}(T-2,\|r\|),\quad
\bar q=\frac{\bar s}{\|\bar s\|}\operatorname{thnt}(T-3,\|\bar s\|).
\]
E114–118은 `bar s=0`이면 `bar q=0`으로 정의한다. 여기에도 같은 영점 처리를 한다.
반면 `r=0`이면 E108의 투영 방향이 정의되지 않는다. 임의 방향을 발명하지 않고 오류로 처리한다.
표준의 짝수 Q 대칭 등간격 격자는 1/2을 포함하지 않지만, 사용자 지정 격자에는 이 전제가
자동 보장되지 않으므로 생성된 codeword 검산에서 확인해야 한다.

### 4.4 E121–134: 평행 성분, contraction, principal square root

\[
v=\begin{cases}
0,&u=0,\\
\dfrac{u}{\|u\|}\sqrt{1-e^{-\|u\|^2}}
 \sqrt{(1-\|p\|^2)(1-\|\bar q\|^2)},&u\ne0,
\end{cases}
\qquad W=[p,\bar q+v],\qquad X=\begin{bmatrix}(I_2-W^HW)^{1/2}\\W\end{bmatrix}.
\]
`v`는 p와 평행, `bar q`는 p에 직교한다. `a=||p||`, `b=||bar q||`, `c=||v||`라 하면
\[
\det(I-W^HW)=(1-a^2)(1-b^2)-c^2
=(1-a^2)(1-b^2)e^{-\|u\|^2}>0
\]
이므로 내부의 유한 입력에서는 위쪽 블록은 positive definite다.
principal Hermitian square root를 A라 하면 `X^HX=A²+W^HW=I₂`다.
Cholesky나 열별 정규화로 A를 바꾸면 참조 행렬이 달라진다. 포트는 eig 재구성을 사용한다.

### 4.5 D54–58, D104–108: subspace 추정과 polar gauge 제거

Y의 상위 두 left singular vector로 `C=[C1;C2]`를 만든다.
`C1=U1 S1 V1^H`, `Q=U1 V1^H`, `W=C2 Q^H`가 D104–108이다.
무잡음에서 `rank(H)=2`이면 C는 `X R` 형태이고 R은 unitary다.
위쪽 A가 positive definite이므로 `polar(A R)=R`; 따라서 `C2 Q^H=W`가 복원된다.
이는 부호 두 개만 맞추는 과정이 아니라 임의의 2×2 unitary basis mixing을 제거하는 과정이다.
G 경로는 같은 subspace를 `eigh(YY^H)`로 얻는다. 충분통계량 주장은 STATE/R1을 인용한다.

### 4.6 D113–130: inverse radial map

\[
p=W_{:,1},\quad t=W_{:,2},\quad v=\frac{p^Ht}{\|p\|^2}p,
\quad \bar q=t-v,
\quad D_0=\sqrt{(1-\|p\|^2)(1-\|\bar q\|^2)},\quad a_0=\|v\|/D_0.
\]
D119–121의 비음수 근과 방향은
\[
\|u\|=\sqrt{-\log(1-a_0^2)},\qquad u=\|u\|\,v/\|v\|.
\]
D124–127의 두 근은 §4.2의 역함수로 얻어
\[
r=\frac{p}{\|p\|}\sqrt{P^{-1}(T-2,I_{\|p\|^2}(T-2,2))},\quad
\bar s=\frac{\bar q}{\|\bar q\|}\sqrt{P^{-1}(T-3,I_{\|\bar q\|^2}(T-3,2))},\quad
s=\bar s+u.
\]
**원본 D121,127에는 영점 가드가 없다.** `v=0`이면 `u=0`, `bar q=0`이면 `bar s=0`이라는
연속극한 처리를 신규 확장으로 명시한다. `u/v`의 크기 비는 v→0에서 `1/D0`이다.
이 확장을 원본에서 이미 구현했다고 적지 않는다. `p=0`, singular C1, 정확히 0인 D0는
다른 비식별 경계이므로 오류로 처리하고, 무잡음 검산에서 그런 블록을 제외하거나 재추출하지 않는다.

### 4.7 D135–163: CDF와 hard decision

복원한 실수 성분 z에 대해 `x_hat=Phi_N(sqrt(2)*z)`를 계산하고
\[
\widehat d_j=\arg\min_{0\le q<Q}(\widehat x_j-\operatorname{lattice}[q])^2
\]
로 결정한다. D143–157의 `min`은 동점에서 첫 MATLAB 인덱스를 고른다.
NumPy `argmin`의 첫 인덱스 규칙을 쓰면 0-기반으로 같다.
`np.rint`는 half-to-even이고, R6의 `floor(Q*x_hat)`는 이 alpha 격자에 일반적으로 맞지 않는다.
따라서 두 방식 모두 가져오지 않는다. D159–163의 순서로 digits를 만들고 R6식 인덱스로 변환한다.

### 4.8 BG25–28, GB25–32, F113–116·128–130: 비트 라벨

각 좌표의 B0자리 자연 이진 정수 q에 대해 Gray 정수는 `g=q XOR (q>>1)`다.
좌표 안에서 MSB가 먼저 나오며, 복원은 Gray bit의 prefix XOR다.
F114–115는 송신 payload가 Gray 라벨이라는 것을 보여 준다:
`payload → Gray2Bin → bit2int → +1 → encoder`.
수신은 `symbol-1 → int2bit → Bin2Gray`다. payload에 Bin2Gray를 먼저 적용하면 방향이 틀린다.

전역 codeword 번호는 R6처럼 `k=sum_j d_j Q^j`로 정한다.
이는 MATLAB 인코더가 정의한 전역 번호가 아니라 **새로운 R6 통합 규약**이다.
따라서 첫 좌표는 k의 최하위 자리이지만, 각 좌표의 label bit는 MSB 우선이다.
`k`의 단순 이진 표현과 `bits[k]`를 같다고 가정해서는 안 된다.

## 5. MATLAB → NumPy에서 반드시 구분할 것

| 지점 | 이식 규약 |
|---|---|
| 1-기반 인덱스 | MATLAB symbol = digit+1. Python 배열 내부·최종 codeword index는 모두 0-기반. +1은 참조 호출 경계에서만 |
| BG26/GB26의 reshape | 실제 shape는 `(B0,nd)`이며 **열 하나가 한 좌표**다. 원본의 “rows” 주석을 따라 뒤집지 않는다 |
| 열 우선 reshape | 한 블록에서 `b.reshape((B0,nd),order='F').T`는 `b.reshape((nd,B0))`와 같다. 배치는 `(n,nd,B0)`의 C-order로 표현한다. `(n,...)` 전체를 order='F'로 reshape해 블록을 섞지 않는다 |
| E/D의 `(:)` | 각 블록 내부의 벡터화다. 배치 축을 포함한 전역 flatten과 다르다. A1에는 A2의 별도 reshape 절차를 이식하지 않는다 |
| `'`와 `.'` | MATLAB `'`는 conjugate transpose. 배치에서는 `.conj().swapaxes(-1,-2)`. 3차원 `.T`는 배치 축까지 뒤집으므로 쓰지 않는다 |
| SVD 반환 | MATLAB `[U,S,V]` 대 NumPy `(U,s,Vh)`. Q는 `U @ Vh`이지 `U @ Vh.conj().T`가 아니다 |
| SVD 크기·순서 | 원본 D57은 full SVD. N>=2이면 reduced SVD의 상위 두 열로 충분하다. N<2에서는 원본 full U의 나머지 열이 식별된 신호 subspace라는 보장이 없어 G1_Y에서 거부한다 |
| eig/eigh 순서 | NumPy eigh는 오름차순이므로 G의 마지막 두 열을 역순으로 취한다. sqrtm용 eig 재구성은 고유값과 고유벡터를 같이 쓰면 정렬 자체에 무관하다 |
| 부호·복소 위상·중복값 | 서로 같은 U 원소를 요구하지 않는다. 선택된 두 차원 안의 임의 unitary 변화는 polar gauge에서 소거된다. 2번째와 3번째 고유값의 경계 동률 또는 singular C1에서는 subspace/gauge가 비유일하며 참조 인덱스 일치를 보장하지 않는다 |
| sqrtm | Hermitian principal square root를 eigenvalue 함수로 재구성한다. Cholesky 대체, 비켤레 전치, elementwise sqrt는 틀린 이식이다 |
| norm과 reduction | MATLAB은 한 블록의 vector norm. NumPy는 각 블록의 좌표 축만 줄인다. 축 없는 norm으로 전체 batch를 합치지 않는다 |
| fzero | 같은 방정식의 유일근을 특수함수로 구한 신규 수치 구현. 3개 scalar solve를 Python 블록 루프로 감싸지 않는다. 반복 경로와 마지막 부동소수 비트 일치는 주장하지 않는다 |
| 난수 | 같은 seed가 MATLAB/NumPy에서 같은 H,Z를 뜻하지 않는다. 교차 언어 대조에는 실제 Y를 공유한다. M=1 R6 대조만 동일 NumPy RNG 호출열을 비교한다 |

수치 정책: double precision에서 unit interval과 PSD를 확인한다. 기계 epsilon의 상수배 이내의
이탈만 clip한다. 역함수의 반올림된 1은 `nextafter(1,0)`으로 제한한다. 이 endpoint 정책도
`NEW`이며 원본의 경계 fzero 동작과 같다고 주장하지 않는다. 큰 위반·비유한 값은 오류다.
원본 D124–125의 [0,100] 밖 radial inverse 역시 오류로 처리한다.

## 6. 코드 — 추출 대상 `code/m2_gate.py`

아래 블록은 한 모듈이다. 원본 R6를 같은 `code/` 경로에서 import한다.
라이선스 고지를 포함해 그대로 추출한다. 라이브러리의 loop는 사전계산 chunk와 bit 축뿐이다.

```python
# NEW: A1 NumPy port; preserve the upstream notice below when extracting this module.
# BSD 3-Clause License
#
# Copyright (c) 2025, Diego Cuevas
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np  # NEW: NumPy batch implementation
from scipy.special import ndtr, ndtri, gammainc, gammaincc, gammaincinv, gammainccinv, betainc, betaincinv, gammaln  # NEW: existing SciPy dependency
import cubesplit_gap as r6  # NEW: unchanged R6 source
from cubesplit_gap import exact_llr, lse, snr_at, gram, LDPC, make_ldpc, gf2_rref, gvec, ml_metrics_fast, exact_llr_chunked  # NEW: reuse, not reimplementation


def _book(C):  # NEW: R6 (T,K) / general (T,M,K) adapter
    C = np.asarray(C)  # NEW
    if C.ndim == 2:  # NEW: preserve the legacy layout
        C = C[:, None, :]  # NEW
    if C.ndim != 3 or not (1 <= C.shape[1] <= C.shape[0]) or C.shape[2] == 0:  # NEW
        raise ValueError('C must have shape (T,K) or (T,M,K), 1<=M<=T, K>0')  # NEW
    if not np.isfinite(C).all():  # NEW
        raise ValueError('nonfinite codebook')  # NEW
    return C  # NEW


def channel(rng, X, N, rho):  # NEW: R6 signature unchanged; spec 3/6 normalization
    X = np.asarray(X)  # NEW
    if X.ndim == 2:  # NEW: legacy rows
        X = X[:, :, None]  # NEW
    if X.ndim != 3 or not (1 <= X.shape[2] <= X.shape[1]):  # NEW
        raise ValueError('X must be (n,T) or (n,T,M), 1<=M<=T')  # NEW
    if not isinstance(N, (int, np.integer)) or N < 1 or np.ndim(rho) != 0 or not np.isfinite(rho) or rho < 0:  # NEW
        raise ValueError('N>=1 integer; rho>=0 finite scalar')  # NEW
    if not np.isfinite(X).all():  # NEW
        raise ValueError('nonfinite X')  # NEW
    n, T, M = X.shape  # NEW: generalizes R6 channel
    H = (rng.standard_normal((n, M, N)) + 1j * rng.standard_normal((n, M, N))) / np.sqrt(2)  # NEW: shape only; R6 draw order
    Z = (rng.standard_normal((n, T, N)) + 1j * rng.standard_normal((n, T, N))) / np.sqrt(2)  # NEW: unchanged R6 noise line
    if M == 1:  # NEW: preserve R6 multiplication order and RNG stream exactly
        return np.sqrt(rho * T) * X[:, :, 0, None] * H[:, None, 0, :] + Z  # NEW: R6 channel
    return np.sqrt(rho * T / M) * (X @ H) + Z  # NEW: spec 3; not a change to the noise


def kappa_of(rho, T, M=1):  # NEW: R6 leading arguments and default M=1
    if not isinstance(T, (int, np.integer)) or not isinstance(M, (int, np.integer)) or not 1 <= M <= T:  # NEW
        raise ValueError('integer dimensions 1<=M<=T')  # NEW
    rho = np.asarray(rho)  # NEW: scalar or SNR array
    if not np.isfinite(rho).all() or np.any(rho < 0):  # NEW
        raise ValueError('rho must be finite and nonnegative')  # NEW
    a = rho * T / M  # NEW: substitute a=rho*T/M in STATE/R1, sigma^2=1
    return a / (1.0 + a)  # NEW: spec 3; M=1 is R6


def codebook_features(C):  # NEW: R6 name, arguments, float32 return retained
    C = _book(C)  # NEW
    T, M, K = C.shape  # NEW
    if M == 1:  # NEW: exact R6 branch, including arithmetic/dtype
        return r6.codebook_features(C[:, 0, :])  # NEW: unchanged implementation reused
    iu = np.triu_indices(T, 1)  # NEW: R6 upper-triangle order
    P = np.einsum('tmk,smk->tsk', C.conj(), C)  # NEW: P is conjugate(XX^H), NOT XX^H
    diag = np.real(np.einsum('ttk->tk', P)).T  # NEW: R6 feature order
    cross = P[iu[0], iu[1], :].T  # NEW: must pair with R6 gvec's -2 Im(G)
    return np.concatenate([diag, cross.real, cross.imag], 1).astype(np.float32)  # NEW: (K,T^2)


def _chordal_rows(C32, start, stop):  # NEW: batched chordal closeness; C32 is (T,M,K)
    T, M, K = C32.shape  # NEW
    left = np.ascontiguousarray(C32[:, :, start:stop].conj().transpose(2, 1, 0))  # NEW
    z = left.reshape(-1, T) @ C32.reshape(T, M * K)  # NEW: (rows*M, M*K)
    return (np.abs(z.reshape(stop - start, M, M, K)) ** 2).sum(axis=(1, 2))  # NEW: ||Xi^H Xj||_F^2


def precompute_lists_sorted(cs, C, bits, eta_max, chunk=256, topm=2048, log=None):  # NEW: exact R6 signature/return structure
    C = _book(C)  # NEW
    bits = np.asarray(bits)  # NEW
    if bits.ndim != 2 or bits.shape[0] != C.shape[2] or not np.isin(bits, (0, 1)).all():  # NEW
        raise ValueError('bits must be (K,B) binary labels')  # NEW
    K, B = bits.shape  # NEW: R6
    if any(not isinstance(v, (int, np.integer)) or v < 1 for v in (eta_max, chunk, topm)):  # NEW: scalar settings loop, not block loop
        raise ValueError('eta_max, chunk, topm must be positive integers')  # NEW
    if eta_max > min(K, topm):  # NEW
        raise ValueError('eta_max must not exceed min(K,topm)')  # NEW
    counts = bits.sum(0)  # NEW
    if np.any(np.minimum(counts, K - counts) < eta_max):  # NEW
        raise ValueError('each bit value needs at least eta_max codewords')  # NEW
    lists = np.empty((K, B, 2, eta_max), np.int32)  # NEW: R6 return dtype/layout
    C32 = np.ascontiguousarray(C, dtype=np.complex64)  # NEW: preserve R6 closeness precision
    mask1 = bits.astype(bool)  # NEW: R6
    nfall = 0  # NEW: number of fallback (center,bit,value) rows, not unique centers
    for s in range(0, K, chunk):  # NEW: chunk loop, never a received-block loop
        corr = _chordal_rows(C32, s, min(s + chunk, K))  # NEW: replaces R6 rank-one closeness
        r = corr.shape[0]  # NEW: R6
        if topm < K:  # NEW: R6 stage 1 retained
            top = np.argpartition(-corr, topm - 1, axis=1)[:, :topm]  # NEW: R6 top pool
            top = np.sort(top, axis=1)  # NEW: canonical index order before stable score sort
            ctop = np.take_along_axis(corr, top, 1)  # NEW: R6
            cut = ctop.min(1)  # NEW: detect ties crossing the topm boundary
            cut_tie = (corr == cut[:, None]).sum(1) > (ctop == cut[:, None]).sum(1)  # NEW
        else:  # NEW: R6
            top = np.broadcast_to(np.arange(K), (r, K))  # NEW: R6 equivalent without tile copy
            ctop = corr  # NEW: R6
            cut = np.full(r, -np.inf)  # NEW
            cut_tie = np.zeros(r, dtype=bool)  # NEW
        btop = mask1[top]  # NEW: R6 stage 2 retained
        for j in range(B):  # NEW: bit-coordinate loop, not received-block loop
            for b in (0, 1):  # NEW: R6 bit-value loop
                v = np.where(btop[:, :, j] == bool(b), ctop, -1.0)  # NEW: R6 mask
                idx = np.argsort(-v, axis=1, kind='stable')[:, :eta_max]  # NEW: R6 stable score sort
                lists[s:s + r, j, b] = np.take_along_axis(top, idx, 1)  # NEW: R6
                last = np.take_along_axis(v, idx, 1)[:, -1]  # NEW: R6 unfilled-list check
                bad = (last < 0) | (cut_tie & (last == cut))  # NEW: also repair truncated ties
                if bad.any():  # NEW: batch fallback replaces R6 loop over rr
                    nfall += int(bad.sum())  # NEW: retain operational fallback-count convention
                    vf = np.where(mask1[:, j] == bool(b), corr[bad], -1.0)  # NEW: full candidate row
                    lists[s:s + r, j, b][bad] = np.argsort(-vf, axis=1, kind='stable')[:, :eta_max]  # NEW: batched stable fallback
        if log is not None and (s // chunk) % 32 == 0:  # NEW: R6 logging cadence
            log(f'    lists {min(s + chunk, K)}/{K}  fallbacks so far {nfall}')  # NEW: R6
    return lists, nfall  # NEW: R6 return unchanged


def _prob(x, open_upper=False):  # NEW: floating-point guard; never repairs material domain violations
    x = np.asarray(x, dtype=np.float64)  # NEW
    tol = 64 * np.finfo(np.float64).eps  # NEW: arithmetic tolerance, not a channel parameter
    if not np.isfinite(x).all() or np.any(x < -tol) or np.any(x > 1 + tol):  # NEW
        raise FloatingPointError('outside the unit interval beyond roundoff')  # NEW
    return np.clip(x, 0.0, np.nextafter(1.0, 0.0) if open_upper else 1.0)  # NEW: documented endpoint policy


def thnt(order, t):  # H1,34-35; NEW: batched evaluation of the same unique root
    t = np.asarray(t, dtype=np.float64)  # NEW
    if not isinstance(order, (int, np.integer)) or order < 1 or not np.isfinite(t).all() or np.any(t < 0):  # NEW
        raise ValueError('positive integer order and finite nonnegative radii required')  # NEW
    p = gammainc(order, t * t)  # NEW(H34): regularized lower gamma is the RHS
    q = gammaincc(order, t * t)  # NEW(H34): complementary tail avoids 1-p cancellation
    z = np.where(p <= 0.5, betaincinv(order, 2, p), 1 - betaincinv(2, order, q))  # NEW(H34-35): z=theta^2
    return np.sqrt(_prob(z))  # NEW(H35): same root, not MATLAB fzero iterates


def _radial(z, order, inverse=False):  # NEW: shared batched radial map and its continuous zero limit
    radius = np.linalg.norm(z, axis=1)  # E107,113; D124-127
    if inverse:  # NEW: algebraic solution of D124-125
        w2 = _prob(radius * radius, open_upper=True)  # NEW(D124-125): roundoff-only guard
        p = betainc(order, 2, w2)  # NEW(H34,D124-125)
        q = betainc(2, order, 1 - w2)  # NEW: stable complementary beta tail
        target = np.sqrt(np.where(p <= 0.5, gammaincinv(order, p), gammainccinv(order, q)))  # NEW(D124-125): same unique roots
        if not np.isfinite(target).all() or np.any(target > 100):  # NEW: retain D124-125's finite root-bracket domain
            raise FloatingPointError('radial inverse outside MATLAB [0,100] domain')  # NEW
        zero_ratio = np.exp(gammaln(order + 2) / (2 * order))  # NEW: inverse continuous limit
    else:  # NEW
        target = thnt(order, radius)  # E112,117
        zero_ratio = np.exp(-gammaln(order + 2) / (2 * order))  # NEW: forward continuous limit
    ratio = np.divide(target, radius, out=np.full_like(radius, zero_ratio), where=radius != 0)  # NEW: avoids eager 0/0 in np.where
    return z * ratio[:, None]  # E112,115,117; D126-127 with NEW zero extension


def _sqrt_psd(A):  # NEW: batched Hermitian principal square root for E134
    A = (A + A.conj().swapaxes(-1, -2)) / 2  # NEW: remove skew-Hermitian roundoff only
    values, vectors = np.linalg.eigh(A)  # NEW(E134): ascending eigenvalues are paired with their vectors
    tol = 64 * np.finfo(np.float64).eps  # NEW: A=I-W^H W has unit scale
    if np.any(values < -tol):  # NEW
        raise FloatingPointError('I-W^H W is not positive semidefinite')  # NEW
    return (vectors * np.sqrt(np.maximum(values, 0))[:, None, :]) @ vectors.conj().swapaxes(-1, -2)  # NEW(E134): principal sqrtm, not Cholesky


class GrassLattice:  # NEW: R6-style A1 interface; no fictitious Cube-split cell
    def __init__(self, T, B0, lattice):  # NEW: B0 is spec B-prime; lattice is required
        if not isinstance(T, (int, np.integer)) or T < 4 or not isinstance(B0, (int, np.integer)) or B0 < 1:  # NEW
            raise ValueError('A1 needs integer T>=4, B0>=1')  # NEW: H requires T-3>=1
        self.T, self.M, self.B0 = int(T), 2, int(B0)  # NEW: M=2 only for this encoder
        self.nd, self.Q = 4 * (self.T - 2), 1 << self.B0  # E53,87-90; F87,91
        self.B = self.nd * self.B0  # F87
        if self.B > 62:  # NEW: integer-index guard, not a claim about feasible exhaustive enumeration
            raise ValueError('index exceeds signed int64 support')  # NEW
        self.K = 1 << self.B  # NEW: R6 indexing interface
        self.grid = np.asarray(lattice, dtype=np.float64).copy()  # E55-57
        if self.grid.shape != (self.Q,) or not np.isfinite(self.grid).all() or np.any(self.grid <= 0) or np.any(self.grid >= 1) or np.any(np.diff(self.grid) <= 0):  # NEW
            raise ValueError('lattice must have Q strictly increasing points in (0,1)')  # NEW
        self.weights = np.left_shift(np.int64(1), self.B0 * np.arange(self.nd, dtype=np.int64))  # NEW: R6 coordinate 0 is least significant

    def _digits(self, digits):  # NEW: validate without silent integer truncation
        digits = np.asarray(digits)  # NEW
        if digits.ndim != 2 or digits.shape[1] != self.nd or digits.dtype.kind not in 'iu' or np.any(digits < 0) or np.any(digits >= self.Q):  # NEW
            raise ValueError('digits must be integer (n,nd) in [0,Q)')  # NEW
        return digits.astype(np.int64, copy=False)  # NEW: Python digit = MATLAB symbol - 1

    def index(self, digits):  # NEW: R6 naming, no cell argument
        return self._digits(digits) @ self.weights  # NEW: k=sum_j digit_j Q^j

    def from_index(self, k):  # NEW: R6 naming
        k = np.asarray(k)  # NEW
        if k.ndim != 1 or k.dtype.kind not in 'iu' or np.any(k < 0) or np.any(k >= self.K):  # NEW
            raise ValueError('indices must be integer (n,) in [0,K)')  # NEW
        return (k.astype(np.int64)[:, None] // self.weights[None, :]) % self.Q  # NEW: batch base-Q expansion

    def bits(self, digits):  # BG25-28; F129-130; NEW: R6 naming
        digits = self._digits(digits)  # NEW
        gray = digits ^ (digits >> 1)  # BG27: adjacent binary-bit XOR, independently for each coordinate
        shifts = self.B0 - 1 - np.arange(self.B0)  # NEW(F129): MSB first within each coordinate
        return ((gray[:, :, None] >> shifts) & 1).astype(np.uint8).reshape(len(digits), self.B)  # BG26,28: per-block column-major groups represented as (n,nd,B0)

    def from_bits(self, bits):  # GB25-32; F114-115; NEW: R6 naming
        bits = np.asarray(bits)  # NEW
        if bits.ndim != 2 or bits.shape[1] != self.B or not np.isin(bits, (0, 1)).all():  # NEW
            raise ValueError('bits must be binary (n,B)')  # NEW
        groups = bits.astype(np.uint8).reshape(len(bits), self.nd, self.B0)  # GB26: per-block MATLAB column groups, not a global Fortran reshape
        binary = np.bitwise_xor.accumulate(groups, axis=2)  # GB27-31: prefix XOR instead of bit-position loop
        shifts = self.B0 - 1 - np.arange(self.B0)  # NEW(F115): bit2int MSB convention
        return (binary.astype(np.int64) << shifts).sum(2)  # F115: omit MATLAB's +1 internally

    def symbols(self, digits):  # E1,53-57; NEW: batch return (n,T,2)
        digits = self._digits(digits)  # NEW
        xtilde = self.grid[digits]  # E55-57: repmat/sub2ind replaced with direct indexing
        d = self.T - 2  # E53
        r = (ndtri(xtilde[:, :2*d:2]) + 1j * ndtri(xtilde[:, 1:2*d:2])) / np.sqrt(2)  # E87-88,96-97,101
        s = (ndtri(xtilde[:, 2*d::2]) + 1j * ndtri(xtilde[:, 2*d+1::2])) / np.sqrt(2)  # E89-90,98-99,102
        nr2 = np.sum(np.abs(r) ** 2, axis=1)  # E107
        if np.any(nr2 == 0):  # NEW: E108 is undefined at r=0; do not invent a projection direction
            raise ValueError('lattice/digits produce r=0, outside the reference chart')  # NEW
        u = ((r.conj() * s).sum(1) / nr2)[:, None] * r  # E108
        sbar = s - u  # E109
        p = _radial(r, self.T - 2)  # E112
        qbar = _radial(sbar, self.T - 3)  # E113-118: includes exact zero branch
        nu = np.linalg.norm(u, axis=1)  # E121,124
        f = np.divide(np.sqrt(-np.expm1(-nu * nu)), nu, out=np.ones_like(nu), where=nu != 0)  # E121-124; NEW: stable expm1 and zero limit
        den = np.sqrt(_prob(1 - np.sum(np.abs(p)**2, 1)) * _prob(1 - np.sum(np.abs(qbar)**2, 1)))  # E125; NEW: roundoff guard
        v = u * (f * den)[:, None]  # E122,125
        W = np.stack((p, qbar + v), axis=2)  # E129
        A = np.eye(2) - W.conj().swapaxes(1, 2) @ W  # E134
        return np.concatenate((_sqrt_psd(A), W), axis=1)  # E134

    def codebook(self):  # NEW: R6 tuple contract, C=(T,M,K), bits=(K,B)
        digits = self.from_index(np.arange(self.K, dtype=np.int64))  # NEW
        C = np.ascontiguousarray(self.symbols(digits).transpose(1, 2, 0))  # NEW: no flatten/reshape of the batch axis
        return C, self.bits(digits)  # NEW: R6-style tuple


def _decode_basis(gl, C):  # NEW: shared batched G1 after D57-58
    left, singular, vh = np.linalg.svd(C[:, :2, :], full_matrices=False)  # D104,106; NumPy returns V^H
    if np.any(singular[:, -1] == 0):  # NEW: polar gauge is nonunique at singular C1
        raise FloatingPointError('singular upper block: reference chart is not unique')  # NEW
    Q = left @ vh  # D107: u*v'; do not conjugate vh again
    W = C[:, 2:, :] @ Q.conj().swapaxes(1, 2)  # D105,108
    p, t = W[:, :, 0], W[:, :, 1]  # D113-114
    p2 = np.sum(np.abs(p)**2, 1)  # D115
    if np.any(p2 == 0):  # NEW: no fabricated direction at p=0
        raise FloatingPointError('p=0, outside the reference inverse chart')  # NEW
    v = ((p.conj() * t).sum(1) / p2)[:, None] * p  # D115
    qbar = t - v  # D116
    den = np.sqrt(_prob(1 - p2) * _prob(1 - np.sum(np.abs(qbar)**2, 1)))  # D120
    if np.any(den == 0):  # NEW: exact singular boundary is not silently assigned an index
        raise FloatingPointError('singular radial inverse denominator')  # NEW
    nv = np.linalg.norm(v, axis=1)  # D120-121
    ratio2 = _prob((nv / den)**2, open_upper=True)  # NEW(D119-120): same scalar equation with roundoff guard
    uu = np.sqrt(-np.log1p(-ratio2))  # NEW(D119-120): analytic nonnegative root
    scale = np.divide(uu, nv, out=1 / den, where=nv != 0)  # NEW(D121): continuous zero extension
    u = v * scale[:, None]  # D121; NEW extension at v=0
    r = _radial(p, gl.T - 2, inverse=True)  # NEW(D124,126): same root, vectorized special function
    sbar = _radial(qbar, gl.T - 3, inverse=True)  # NEW(D125,127): includes qbar=0 limit
    s = sbar + u  # D130
    coords = np.concatenate((np.stack((r.real, r.imag), 2).reshape(len(C), -1), np.stack((s.real, s.imag), 2).reshape(len(C), -1)), 1)  # D135-138,159-163: preserve coordinate order
    estimates = ndtr(np.sqrt(2) * coords)  # D135-138: N(0,1/2) CDF
    digits = np.argmin((estimates[:, :, None] - gl.grid[None, None, :])**2, axis=2)  # D143-163: nearest level; first index wins exact ties
    return gl.index(digits)  # NEW: zero-based (n,) codeword index, not MATLAB coordinate vector


def decode_G1_Y(gl, Y):  # NEW: explicit Y API, avoids N=T input ambiguity
    Y = np.asarray(Y, dtype=np.complex128)  # NEW: match MATLAB double precision
    if Y.ndim != 3 or Y.shape[1] != gl.T or Y.shape[2] < 2 or not np.isfinite(Y).all():  # NEW
        raise ValueError('Y must be finite (n,T,N), N>=2')  # NEW
    U, singular, _ = np.linalg.svd(Y, full_matrices=False)  # D54,57: top two columns suffice when N>=2
    if np.any(singular[:, 1] <= 64 * np.finfo(float).eps * singular[:, 0]):  # NEW: reject numerically rank-deficient inputs
        raise FloatingPointError('received matrix has numerical rank below 2')  # NEW
    return _decode_basis(gl, U[:, :, :2])  # D58; NEW: batched index output


def decode_G1(gl, G):  # NEW: R6-like (codebook,G) API; output is (n,) indices
    G = np.asarray(G, dtype=np.complex128)  # NEW
    if G.ndim != 3 or G.shape[1:] != (gl.T, gl.T) or not np.isfinite(G).all():  # NEW
        raise ValueError('G must be finite (n,T,T)')  # NEW
    tol = 64 * np.finfo(float).eps * np.maximum(np.linalg.norm(G, axis=(1, 2)), np.finfo(float).tiny)  # NEW: scale-aware arithmetic tolerance
    if np.any(np.linalg.norm(G - G.conj().swapaxes(1, 2), axis=(1, 2)) > tol):  # NEW
        raise ValueError('G must be Hermitian')  # NEW
    values, vectors = np.linalg.eigh(G)  # NEW(D57-58): G=YY^H has the same principal left subspace
    if np.any(values[:, 0] < -tol) or np.any(values[:, -2] <= tol):  # NEW: PSD and numerical-rank checks
        raise FloatingPointError('G must be PSD with numerical rank at least 2')  # NEW
    return _decode_basis(gl, vectors[:, :, -2:][:, :, ::-1])  # NEW(D58): eigh ascending, take largest two in descending order
```

## 7. 검산 방법

이 절은 **실행할 assertion과 범위**다. 성능 추정치나 벤치마크를 쓰지 않는다.

| EXP §7 | 이 문서로 확인하는 것 / 남는 것 |
|---|---|
| 1 | R6 Cube-split 입력에 대해 `(n,T)`와 `(n,T,1)`를 각각 비교. Y, RNG 잔여 상태, κ, Phi, 지표, ML 결정, 오류 개수 및 exact LLR까지 대조. 지표 상대오차 <1e-6과 동일 SER를 assert |
| 2 | B0=1,2의 모든 codeword에서 `||X^HX-I₂||F<1e-10`. 샘플만 보는 것이 아님 |
| 3 | 1000 무작위 인덱스 왕복. 별도의 bit-group XOR 구성으로 라벨 순서도 검사하여 잘못된 encoder/decoder가 서로 상쇄되는 경우를 막음 |
| 4 | G1에 대해 2000 무잡음 블록을 Y 경로와 G 경로로 검사. **둘 다 모든 인덱스가 일치해야 함**. G2 부분은 다음 chat |
| 5 | G2 미구현이므로 해당 없음. G1 점수 검사를 G2 검산이라고 바꾸어 부르지 않음 |
| 6 | 이 환경에 Octave가 없어 미실행. 아래 참조 대조 절차로 따로 확인해야 함 |
| 7 | 독립 complex128 trace 지표와 기존 float32 gvec/Phi 지표의 argmax를 대조. float32 Phi를 float64로 cast한 값을 ground truth로 쓰지 않음. 일치율 기준은 사양의 ≥0.999 |
| 8 | 같은 float32 closeness의 eta번째와 eta+1번째 경계가 정확히 같은 횟수/분모를 기록. B=8이 기본이고 B=16 전체 setup은 명시적 full_lists 옵션으로 실행. skip을 전체 통과라고 쓰지 않음 |
| 9 | LDPC는 재작성하지 않음. 원본 `sanity.py`와 R6 설정의 구성 검사를 별도로 실행. 여기서 새로 인증하지 않음 |

### 7.1 특히 item 1: M=1 환원

송신 인덱스 생성용 RNG와 채널 RNG를 분리한다. 동일 seed로 만든 독립 RNG 세 개에 대해
기존 R6 channel, 신규 channel의 legacy 입력, 신규 channel의 M=1 입력을 각각 호출한다.
같은 출력만이 아니라 `bit_generator.state`도 비교해 몰래 추가된 난수 소모를 잡는다.
R6 codebook, R6 gvec, R6 ML, R6 exact LLR을 그대로 사용한다.
SER는 새로운 수치를 보고하지 않고 동일 송신 k에 대한 **오류 개수의 equality**로 검산한다.
small Cube-split fixture의 통과는 모든 R6 실행 설정을 재현했다는 주장이 아니다.

### 7.2 특히 item 4: 무잡음 복원

각 B0에서 2000개 k를 고정 seed로 뽑고 `X=gl.symbols(gl.from_index(k))`를 만든다.
N>=2의 복소 Gaussian H를 한 번 생성하고 rank를 확인한다.
**Y=sqrt(rho*T/2)*(X@H), Z=0**를 직접 만든다. `channel`은 항상 잡음을 더하므로 이 검산에 호출하지 않는다.
`decode_G1_Y(gl,Y)==k`와 `decode_G1(gl,gram(Y))==k`를 모두 assert한다.
하나라도 틀리거나 예외가 나면 검산 실패이며, 실패 블록을 버리거나 H를 다시 뽑아 분모를 바꾸지 않는다.
아울러 전체 codebook을 고정 full-rank H로 복원하여 평행/직교 성분이 0인 경우도 검사한다.
이것은 물리적 `N=M` 특수 경우를 주장하는 것이 아니라 대수적인 fixture다.

### 7.3 추출 대상 `code/sanity_m2.py`

다음 코드는 새 검산 코드다. production 함수는 모두 앞 블록이나 원본 R6를 호출한다.
JSON은 YAML의 부분집합이므로 설정 파일은 JSON 표기의 YAML로 저장하고 표준 라이브러리로 읽는다.
별도 YAML 의존성을 설치하거나 기존 설정 loader가 있다고 가정하지 않는다.

```python
import json  # NEW: JSON is the supported dependency-free subset of YAML
import sys  # NEW
from pathlib import Path  # NEW
import numpy as np  # NEW
import cubesplit_gap as r6  # NEW: import the actual repository module on the server
import m2_gate as m2  # NEW: first Python block of GPT01


def check_m1(cfg):  # NEW: spec 7 item 1; unchanged Cube-split encoder
    rng = np.random.default_rng(cfg['seed'])  # NEW
    cs = r6.CubeSplit(4, 1)  # NEW: small R6 fixture, not a SER sweep
    C, labels = cs.codebook()  # NEW
    k = rng.integers(cs.K, size=2000)  # NEW: fixed input indices, separate from channel RNGs
    X = C[:, k].T  # NEW: legacy (n,T)
    phi0 = r6.codebook_features(C)  # NEW
    for N in cfg['N_m1']:  # NEW: parameter loop, not a block loop
        for rho in cfg['rho_linear']:  # NEW: parameter loop, no performance interpolation
            a, b, c = (np.random.default_rng(cfg['seed']) for _ in range(3))  # NEW: independent identical RNG states
            Y0 = r6.channel(a, X, N, rho)  # NEW
            Y1 = m2.channel(b, X, N, rho)  # NEW
            Y2 = m2.channel(c, X[:, :, None], N, rho)  # NEW: explicit M=1 representation
            assert np.array_equal(Y0, Y1) and np.array_equal(Y0, Y2)  # NEW: stronger than just SER equality
            assert a.bit_generator.state == b.bit_generator.state == c.bit_generator.state  # NEW: detect hidden draws
            assert m2.kappa_of(rho, cs.T) == r6.kappa_of(rho, cs.T)  # NEW
            assert m2.kappa_of(rho, cs.T, 1) == r6.kappa_of(rho, cs.T)  # NEW
            assert np.array_equal(m2.codebook_features(C), phi0)  # NEW
            phi1 = m2.codebook_features(C[:, None, :])  # NEW
            assert np.array_equal(phi1, phi0)  # NEW
            G0, G2 = r6.gram(Y0), r6.gram(Y2)  # NEW: existing gram, no replacement
            met0, met1 = r6.ml_metrics_fast(G0, phi0), r6.ml_metrics_fast(G2, phi1)  # NEW
            rel = np.linalg.norm((met0 - met1).astype(float)) / max(np.linalg.norm(met0.astype(float)), np.finfo(float).tiny)  # NEW
            assert rel < 1e-6  # NEW: spec item 1
            pred0, pred1 = met0.argmax(1), met1.argmax(1)  # NEW
            assert np.array_equal(pred0, pred1)  # NEW: includes a consistent first-index ML tie rule
            assert np.count_nonzero(pred0 != k) == np.count_nonzero(pred1 != k)  # NEW: exactly equal SER, no rates reported
            L0 = r6.exact_llr_chunked(G0, phi0, labels, r6.kappa_of(rho, cs.T))  # NEW: unchanged exact LLR
            L1 = r6.exact_llr_chunked(G2, phi1, labels, m2.kappa_of(rho, cs.T, 1))  # NEW
            assert np.array_equal(L0, L1)  # NEW: additional integration check


def check_a1(cfg):  # NEW: spec 7 items 2,3,4(G1),7; no G2 in this document
    rng = np.random.default_rng(cfg['seed'])  # NEW
    for B0 in (1, 2):  # NEW: spec T=4, B=8 and B=16
        Q = 1 << B0  # NEW
        alpha = cfg['alpha_fixture']  # NEW: explicit fixture, NOT a chosen gate alpha
        lattice = alpha + np.arange(Q) * (1 - 2 * alpha) / (Q - 1)  # F91,99-100
        gl = m2.GrassLattice(4, B0, lattice)  # NEW
        C, labels = gl.codebook()  # NEW
        Xall = C.transpose(2, 0, 1)  # NEW
        error = np.linalg.norm(Xall.conj().swapaxes(1, 2) @ Xall - np.eye(2), axis=(1, 2))  # NEW
        assert np.all(error < 1e-10)  # NEW: item 2, all K, not just transmitted words
        k = rng.integers(gl.K, size=1000)  # NEW: item 3
        d = gl.from_index(k)  # NEW
        assert np.array_equal(gl.index(gl.from_bits(gl.bits(d))), k)  # NEW
        assert np.array_equal(gl.from_bits(labels[k]), d)  # NEW
        binary = ((d[:, :, None] >> (B0 - 1 - np.arange(B0))) & 1).astype(np.uint8)  # NEW: independent bit-group fixture
        gray = np.concatenate((binary[:, :, :1], binary[:, :, 1:] ^ binary[:, :, :-1]), 2)  # BG27; NEW: independent oracle, not a mapper in the library
        assert np.array_equal(gray.reshape(1000, gl.B), labels[k])  # NEW: round-trip alone would miss a shared permutation error
        k0 = rng.integers(gl.K, size=2000)  # NEW: item 4; never drop failures or singular blocks
        X0 = gl.symbols(gl.from_index(k0))  # NEW
        for N in cfg['N_a1']:  # NEW: require N>=2; N=M is not given a special physical interpretation
            assert N >= 2  # NEW
            H = (rng.standard_normal((2000, 2, N)) + 1j * rng.standard_normal((2000, 2, N))) / np.sqrt(2)  # NEW: same law as spec 3
            assert np.all(np.linalg.matrix_rank(H) == 2)  # NEW: do not redraw/discard rank failures
            Y0 = np.sqrt(gl.T / 2) * (X0 @ H)  # NEW: rho=1, Z=0; do NOT call channel(), which always adds noise
            assert np.array_equal(m2.decode_G1_Y(gl, Y0), k0)  # NEW: item 4 G1 in Y form, 2000/2000 required
            assert np.array_equal(m2.decode_G1(gl, r6.gram(Y0)), k0)  # NEW: item 4 G1 in G form
        fixed = np.broadcast_to(np.eye(2), (gl.K, 2, 2))  # NEW: additional all-codeword identity-channel fixture
        assert np.array_equal(m2.decode_G1_Y(gl, Xall @ fixed), np.arange(gl.K))  # NEW: specifically exercises zero parallel/orthogonal components
        phi = m2.codebook_features(C)  # NEW
        ncheck = cfg['precision_blocks']  # NEW: choose before execution; do not change after seeing disagreements
        k1 = rng.integers(gl.K, size=ncheck)  # NEW
        rho = cfg['precision_rho_linear']  # NEW: linear SNR input, not dB
        Y = m2.channel(rng, gl.symbols(gl.from_index(k1)), cfg['N_a1'][-1], rho)  # NEW
        G = r6.gram(Y)  # NEW
        pred32, pred64 = np.empty(ncheck, np.int64), np.empty(ncheck, np.int64)  # NEW
        P = np.einsum('tmk,smk->tsk', C, C.conj())  # NEW: actual XX^H, independent of the feature convention
        for s in range(0, ncheck, cfg['chunk']):  # NEW: bounded batch chunk, not per-block Python loop
            e = min(s + cfg['chunk'], ncheck)  # NEW
            m64 = np.einsum('nts,stk->nk', G[s:e], P).real  # NEW: direct complex128 trace, not float32 Phi cast to float64
            m32 = r6.ml_metrics_fast(G[s:e], phi)  # NEW
            assert np.linalg.norm(m64 - m32) / np.linalg.norm(m64) < 1e-6  # NEW: feature/gvec identity integration
            pred32[s:e], pred64[s:e] = m32.argmax(1), m64.argmax(1)  # NEW
        assert np.mean(pred32 == pred64) >= 0.999  # NEW: item 7
        if cfg['check_lists'] and (B0 == 1 or cfg['full_lists']):  # NEW: B=16 full setup is opt-in, never claim skipped coverage
            etas = cfg['etas']  # NEW
            em = max(etas) + 1  # NEW: need eta+1 to test the boundary
            lists, _ = m2.precompute_lists_sorted(gl, C, labels, em, cfg['chunk'], cfg['topm'])  # NEW: unchanged argument order
            C32 = np.ascontiguousarray(C, dtype=np.complex64)  # NEW
            ties = {eta: 0 for eta in etas}  # NEW: exact computed-score ties only
            for s in range(0, gl.K, cfg['chunk']):  # NEW: batch chunk
                e = min(s + cfg['chunk'], gl.K)  # NEW
                corr = m2._chordal_rows(C32, s, e)  # NEW
                for j in range(gl.B):  # NEW: label-coordinate loop
                    for b in (0, 1):  # NEW
                        scores = np.where(labels[:, j] == b, corr, -1.0)  # NEW: full-row oracle
                        oracle = np.argsort(-scores, axis=1, kind='stable')[:, :em]  # NEW: full stable order only in verification
                        assert np.array_equal(lists[s:e, j, b], oracle)  # NEW: tests top-pool truncation/tie repair
                        selected = np.take_along_axis(corr, oracle, 1)  # NEW
                        for eta in etas:  # NEW: list-size loop
                            ties[eta] += int(np.count_nonzero(selected[:, eta - 1] == selected[:, eta]))  # NEW: item 8
            for eta in etas:  # NEW
                print({'B': gl.B, 'eta': eta, 'boundary_equal_f32': ties[eta], 'denominator': gl.K * gl.B * 2})  # NEW: observed counts only when actually run


def check_generic(cfg):  # NEW: general-M arithmetic checks, not a channel-statistics experiment
    rng = np.random.default_rng(cfg['seed'])  # NEW
    for T in (1, 4):  # NEW: include M=T and the scalar edge case
        for M in range(1, T + 1):  # NEW: rank-parameter loop, not a block loop
            A = rng.standard_normal((16, T, M)) + 1j * rng.standard_normal((16, T, M))  # NEW
            X, _ = np.linalg.qr(A, mode='reduced')  # NEW: arbitrary semiunitary test frames
            C = X.transpose(1, 2, 0)  # NEW: general codebook layout
            Y = m2.channel(rng, X, 2, 1.0)  # NEW: fixture rho=1
            G = r6.gram(Y)  # NEW: unchanged gram
            direct = np.einsum('tmk,nts,smk->nk', C.conj(), G, C).real  # NEW: independent trace
            fast = r6.ml_metrics_fast(G, m2.codebook_features(C))  # NEW
            assert np.linalg.norm(fast - direct) / np.linalg.norm(direct) < 1e-6  # NEW
            labels = ((np.arange(16)[:, None] >> np.arange(3, -1, -1)) & 1).astype(np.uint8)  # NEW: balanced test labels
            corr = m2._chordal_rows(np.ascontiguousarray(C, dtype=np.complex64), 0, 16)  # NEW
            lists, _ = m2.precompute_lists_sorted(None, C, labels, 4, chunk=4, topm=8)  # NEW: force the two-stage/fallback path
            for j in range(4):  # NEW: bit-coordinate loop
                for b in (0, 1):  # NEW
                    oracle = np.argsort(-np.where(labels[:, j] == b, corr, -1.0), axis=1, kind='stable')[:, :4]  # NEW
                    assert np.array_equal(lists[:, j, b], oracle)  # NEW: includes all-equal-subspace case M=T


if __name__ == '__main__':  # NEW: repository-root invocation below
    cfg = json.loads(Path(sys.argv[1]).read_text(encoding='utf-8'))  # NEW: config file is JSON-compatible YAML
    check_m1(cfg)  # NEW
    check_generic(cfg)  # NEW
    check_a1(cfg)  # NEW
    print('requested sanity assertions completed; G2, Octave, LDPC are not certified here')  # NEW: no performance estimate
```

### 7.4 설정과 실행 명령

다음은 **전부 신규 검산 설정**이며 gate 설정이 아니다. `alpha_fixture=0.01`은 F58에서 가져온
검산 입력이지 최적 alpha가 아니다. rho 값들은 정규화 산술을 확인하는 입력이며 측정 결과가 아니다.
`configs/GPT01_sanity.yaml`에 아래 JSON-표기 YAML을 저장한다.

```json
{
  "seed": 11,
  "N_m1": [1, 2, 4],
  "N_a1": [2, 4],
  "rho_linear": [0.0, 1.0, 10.0],
  "alpha_fixture": 0.01,
  "precision_blocks": 1000,
  "precision_rho_linear": 1.0,
  "chunk": 32,
  "etas": [8, 16],
  "topm": 2048,
  "check_lists": true,
  "full_lists": false
}
```

저장소 루트에서 실행한다. 이 문서의 첫 Python 블록을 `code/m2_gate.py`, 둘째 블록을
`code/sanity_m2.py`로 저장한 뒤, 서버의 절대경로 Python을 사용한다.
아래 셸 줄은 모두 신규 실행 절차이며 MATLAB 대응 줄은 없다.

```bash
PY="$HOME/miniforge3/envs/torch/bin/python"  # NEW: project interpreter rule
mkdir -p logs  # NEW: local logs, not committed by this document
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$PY" code/sanity_m2.py configs/GPT01_sanity.yaml > logs/GPT01_sanity.log 2>&1  # NEW: assertion failure preserves exit status
(cd code && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$PY" sanity.py) > logs/GPT01_R6_sanity.log 2>&1  # NEW: unchanged R6/LDPC checks
```

**B=16 전체 리스트를 확인하려면 같은 설정 파일의 `full_lists`를 `true`로 바꾸어 별도 실행한다.**
B=16의 전체 리스트를 실행하지 않은 로그에는 그 검사를 미실행으로 표시한다.
실제 gate lattice가 정해지면 그 입력으로 다시 실행해야 하며 fixture alpha를 gate로 승격하지 않는다.
이 문서는 숫자 로그를 `02_RESULTS_LOG.md`로 옮기거나 실험 완료를 선언하지 않는다.

### 7.5 §8 참조 대조를 나중에 실행할 때

먼저 grassbox의 실제 Git commit을 고정하고, 해당 commit의 E/D/H/BG/GB 파일을 첨부 파일과
비교한다. 아래 첨부 SHA-256을 grassbox Git commit이라고 쓰지 않는다.
Python에서 200개 digits와 실제 Y를 저장하고, 참조에 들어가는 symbol은 `digits+1`로 한다.
MATLAB/Octave의 Y는 각 블록의 `(T,N)` 행렬이며 batch 축을 명시적으로 옮긴다.
원본 인코더 출력은 `(T,2,200)`에서 `(200,T,2)`로 축만 옮겨 비교한다.
원본 복호기의 symbol에서 1을 뺀 후 `sum_j digit_j Q^j`로 Python과 같은 전역 인덱스를 만든다.
codeword 허용오차는 사양의 1e-10, 복호 인덱스는 완전 일치다.

동일 seed의 cross-language RNG에 의존하지 말고 동일 Y 배열을 양쪽에 공급해야 한다.
출력 `.npz`에는 digits, lattice, Y, X_ref, k_ref, 고정한 grassbox commit을 함께 저장한다.
원본이 NaN/실패를 내는 영점 사례는 숨기거나 재추출하지 말고 원본 한계와 NEW 확장의 차이로 기록한다.
단순 NumPy 왕복이나 H34 잔차 검사는 이 참조 대조를 대체하지 않는다.

## 8. ⚠️ 한계와 인수인계

**참조 구현과 대조하지 못했다.** 이 환경에는 Octave가 없었고 MATLAB도 실행하지 않았다.
따라서 “원논문 복호기와 같다”거나 “§7 전체 통과”라고 쓰지 않는다.
여기서 얻은 것은 .m 실행부에서 유도한 수학적으로 대응하는 포트와 명시적 경계 확장이다.
특수함수 근, 영점 연속극한, unit-interval/PSD 수치 가드, numerical-rank 검사,
canonical 리스트 동률 처리는 원본과 구별한 신규 구현이다.
선택 subspace 경계 동률이나 singular chart에서 양 라이브러리의 인덱스 일치를 보장하지 않는다.

이 환경에서는 문서의 Python 코드와 검산 코드를 임시 파일에서 실행해 assertion을 확인했다.
R6 대조에는 연결 도구로 읽은 관련 원본 정의를 옮긴 **로컬 fixture**를 사용했다.
이 fixture는 배포하지 않으며 원본 R6 파일 전체를 import한 서버 실행을 대신하지 않는다.
검산에는 M=1 경로·특징 내적·일반 rank의 무작위 semiunitary 입력·A1의 두 크기·Gray 왕복·
무잡음 G1·float32 대 float64 지표·작은 codebook의 리스트가 포함되었다.
B=16 전체 리스트, Octave 대조, G2, 서버의 LDPC 검산은 이 실행으로 인증하지 않았다.
검산 수치·시간·성능 결과는 여기서 새로 인용하지 않는다.

Claude 인수인계 순서는 코드 두 블록 추출 → 실제 R6 전체 파일을 사용한 검산 →
격자 고정 → 해당 격자로 재검산 → 가능한 환경에서 §8 대조다.
A2·G2 및 G2 중심 soft 비교는 별도 작업이다. 이 제안을 STATE “확정”에 자동으로 올리지 않는다.

## 부록 A. 첨부 원본 고정 정보

아래 값은 이 대화에 실제 첨부된 바이트의 SHA-256이다. upstream Git commit의 증거가 아니다.

| 파일 | SHA-256 |
|---|---|
| GrassLatticeEncoding.m | `71567a7e569b4d692fe29d629a53c36f257a0ac5ccc1d502b0b570262e46124c` |
| GrassLatticeDecoding.m | `8df1e8b27cb3f088f7904c246f2a792f54f60bcc6652cc1c5b994a1954241d59` |
| thnt.m | `8ef5534944b14cd24a768e68cac86d3857509d080a48f47b41254d51a8badd10` |
| Bin2Gray.m | `51f7f718805f55201e9dba97c6f7b9d67fff093c4c0b7935d431a7db580dda13` |
| Gray2Bin.m | `8aa1088898d008a443cdf1d24861001af3de92abf748afe282567e55f86ab249` |
| FindAlphaOpt.m | `44cc9d49813b59bbcbcc5ee3c6af6f4096b2d1dd6e96088066d3f6acdf823688` |
| grassbox_LICENSE.txt | `6a725b2d47dec98c34a98a76239f2277acedd59e29daa45650d5c353de562ba4` |

라이브러리 코드 머리말에 첨부 BSD-3 고지를 보존했다. 이를 코드로 추출할 때도 삭제하지 않는다.
