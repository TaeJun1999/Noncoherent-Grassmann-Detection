# GPT01 — 공통 층 M 일반화와 A1 Grass-Lattice 이식

날짜: 2026-09-18. 상태: **[VERIFY] 제안, 독립 검산 전**.
입력: `EXP_M2_gate_spec.md` §3·4·6·7·8, `cubesplit_gap.py` R6의 해당 함수,
`00_STATE.md` “확정”, `02_RESULTS_LOG.md` R1·R6의 규약, `03_LIT_LOG.md` L1 A1·E.
실제 첨부 사양 파일명은 `EXP_M2_gate_spec (1).md`다.
MATLAB 입력은 `GrassLatticeEncoding.m`, `GrassLatticeDecoding.m`, `thnt.m`,
`Bin2Gray.m`, `Gray2Bin.m`, `FindAlphaOpt.m`의 아래 줄 범위다.
저장소 기준 커밋은 `5fd114a21750354b370525e51d8fe5170e933c35`다.
원본 R6 blob은 `ccf086b1c09e7b393187014817c96f590b222a60`이며 첨부 파일과 일치한다.

## 결론

[VERIFY] 공통 네 함수의 M 일반화와 배치 A1 인코더·Gray 라벨·G1 코드 초안을 제시한다.
특징은 T²를 유지하며, κ는 `(rho*T/M)/(1+rho*T/M)`이다. R6의 gvec과 기존 LLR을 그대로 쓴다.
격자는 필수 입력이다. MATLAB 대조와 실제 gate 격자 검산 전에는 참조 구현과의 동일성을 확정하지 않는다.

## 유도

### 1. 출처, 범위, 코드 계약

사양은 [EXP §3·4·6](../handoff/EXP_M2_gate_spec.md)이다.
posterior와 충분통계량은 [STATE “확정”](../00_STATE.md), R1을 인용만 한다.
문헌 인용은 [LIT L1 A1·E](../03_LIT_LOG.md)로 제한한다.
아래 사상의 직접 근거는 첨부 MATLAB 실행문이다. 논문 본문을 새로 확인했다는 뜻이 아니다.
A2, G2, soft demapper, 성능 스윕과 성능 추정은 이번 산출물에 없다.

줄 번호는 첨부 원본의 빈 줄·주석을 포함한 1-기반이다.
`E=GrassLatticeEncoding.m`, `D=GrassLatticeDecoding.m`, `H=thnt.m`,
`BG=Bin2Gray.m`, `GB=Gray2Bin.m`, `F=FindAlphaOpt.m`으로 표시한다.
`NEW`는 MATLAB에 대응 실행문이 없는 어댑터·검사·일반화다.
`NEW(D124)`는 그 방정식의 **신규 수치 구현**이며 fzero의 직역이 아니다.
모든 Python 실행문에 원본 줄 또는 `NEW`를 붙였다.

| 함수 또는 표현 | 계약 |
|---|---|
| `channel(rng,X,N,rho)` | 기존 `(n,T)` 또는 신규 `(n,T,M)` → `(n,T,N)` |
| `kappa_of(rho,T,M=1)` | 기존 선행 인자 순서 유지. rho는 선형 SNR |
| `codebook_features(C)` | `(T,K)` 또는 `(T,M,K)` → `(K,T²)` float32 |
| `precompute_lists_sorted(cs,C,bits,eta_max,chunk=256,topm=2048,log=None)` | 기존 인자 순서 유지. `(lists,nfall)` 반환 |
| `lists` | `(K,B,2,eta_max)` int32. 자기 자신도 해당 bit-value 집합에 포함 |
| `GrassLattice(T,B0,lattice)` | B0는 사양의 B′. M=2, T≥4 |
| A1 `from_index`, `index` | 각각 `(n,)` → `(n,nd)`, `(n,nd)` → `(n,)` |
| A1 `bits`, `from_bits` | digits와 `(n,B)` uint8 Gray payload 사이의 변환 |
| A1 `symbols`, `codebook` | `(n,T,2)` 및 `(C,bits)` 반환. C의 codeword 축은 마지막 |
| `decode_G1(gl,G)`, `decode_G1_Y(gl,Y)` | 각각 `(n,T,T)`, `(n,T,N)` → 0-기반 `(n,)` int64 |

A1에는 Cube-split cell이 없다. 가짜 cell 인자를 추가하지 않는다.
`N=T`일 때 모양만으로 Y와 G를 구별할 수 없어 입력 종류별 API를 명시했다.
G1에는 N≥2와 수치적으로 식별 가능한 rank-two subspace를 요구한다.
공통 채널은 N≥1을 받는다. 복호 배치는 n≥1을 대상으로 한다.

`exact_llr`, `lse`, `snr_at`, `gram`, `LDPC`, `make_ldpc`, `gf2_rref`,
`gvec`, `ml_metrics_fast`, `exact_llr_chunked`는 원본 import만 한다.
R6의 `metric_of`, `ml_metrics`, `soft_P26_multi`까지 M 일반화됐다는 주장은 하지 않는다.
production 모듈은 `code/m2_gate.py`로 추출하도록 작성했다.

**격자 미지정:** EXP §4는 alpha를 고정하지 않았다.
F91·99–100의 격자는 `lattice[j]=alpha+j*(1-2*alpha)/(Q-1)`이다.
따라서 생성자의 lattice를 필수로 받는다. 임의의 최적 alpha나 midpoint 기본값을 만들지 않는다.
검산의 alpha는 F58에 있는 **입력 fixture**일 뿐, gate 설정이 아니다.
E/D의 lattice 크기 주석과 실행부가 다르다. E55–57·F91·100의 실행부대로 길이는 Q다.

### 2. 공통 층의 수학 — [VERIFY V1]

#### 2.1 채널과 κ의 정규화

사양 §3에서 `a=rho*T/M`라 놓는다.

$$
Y=\sqrt a\,XH+Z,\qquad X^HX=I_M,\qquad a=\frac{\rho T}{M}.
\tag{1}
$$

H와 Z의 원소 분산은 모두 1이다. 따라서

$$
\mathbb E\|\sqrt aXH\|_F^2
=a\,\mathbb E\operatorname{tr}(H^HX^HXH)
=a\,\mathbb E\|H\|_F^2=aMN=\rho TN.
\tag{2}
$$

수신 안테나당 블록 신호 에너지는 rho*T다.
채널 사용당·수신 안테나당 신호 에너지와 잡음 분산의 비가 rho다.
STATE/R1의 확정 식 `κ=γ/[σ²(γ+σ²)]`에 **γ=a, σ²=1**만 대입한다.

$$
\boxed{\kappa=\frac{a}{1+a}=\frac{\rho T/M}{1+\rho T/M}}.
\tag{3}
$$

이는 posterior 재유도가 아니다. M=1이면 a=rho*T이며 R6의 κ가 된다.
또한 `XH=x hᵀ`이므로 신호 생성식도 R6로 환원된다.
코드는 M=1의 난수 호출 순서와 원소별 곱셈 순서를 보존한다.
`10log10(a)=10log10(rho*T)-10log10(M)`는 계수 변환이다.
이 식에서 서로 다른 codebook의 SER 곡선이 일정하게 이동한다고 추론하지 않는다.

#### 2.2 trace와 기존 gvec의 정확한 짝

`P_k=X_kX_k^H`, `⟨A,B⟩=Re tr(A^HB)`로 정의한다.
G가 Hermitian이면 trace 순환성과 실수성에서

$$
\operatorname{tr}(X_k^HGX_k)=\operatorname{tr}(GP_k)=\langle G,P_k\rangle.
\tag{4}
$$

항등식 자체에는 X의 semiunitarity도 필요하지 않다.
다만 R6 특징의 내부 배열은 P_k가 아니라 다음 배열이다.

$$
\widetilde P_{ij,k}=\sum_m\overline{X_{im,k}}X_{jm,k}
=\overline{(P_k)_{ij}}=(P_k)_{ji}.
\tag{5}
$$

R6 순서 그대로 diagonal, upper-real, upper-imag를 연결한다.

$$
\begin{aligned}
\Phi_k&=(\widetilde P_{ii,k};\Re\widetilde P_{ij,k};\Im\widetilde P_{ij,k})_{i<j},\\
g(G)&=(G_{ii};2\Re G_{ij};-2\Im G_{ij})_{i<j}.
\end{aligned}
\tag{6}
$$

상삼각 항과 하삼각 항을 쌍으로 묶으면

$$
\begin{aligned}
\Phi_kg(G)^T
&=\sum_iG_{ii}P_{ii,k}+2\Re\sum_{i<j}G_{ij}\widetilde P_{ij,k}\\
&=\sum_{ij}G_{ij}P_{ji,k}=\operatorname{tr}(GP_k).
\end{aligned}
\tag{7}
$$

특징 차원은 `T+2*T(T-1)/2=T²`다.
**실제 P_k의 상삼각 허수부를 그대로 넣으면 기존 gvec과 부호가 어긋난다.**
코드는 `C.conj(),C` 순서의 einsum을 유지한다.
`np.triu_indices(T,1)`의 순서도 유지하므로 기존 fast ML과 chunked LLR을 재사용한다.

#### 2.3 chordal 근접도와 두 단계 리스트

같은 rank M의 semiunitary codeword에 대해

$$
\begin{aligned}
s_{ij}&=\|X_i^HX_j\|_F^2
=\operatorname{tr}(X_j^HX_iX_i^HX_j)=\operatorname{tr}(P_iP_j),\\
\|P_i-P_j\|_F^2&=\operatorname{tr}(P_i)+\operatorname{tr}(P_j)-2\operatorname{tr}(P_iP_j)
=2M-2s_{ij},\\
d_c^2(X_i,X_j)&=M-s_{ij}=\tfrac12\|P_i-P_j\|_F^2.
\end{aligned}
\tag{8}
$$

따라서 s를 내림차순으로 정렬한다. M=1이면 `s=|c_i^Hc_j|²`다.
정규화 거리 `d_c²/M`도 같은 순서를 준다. 구현에 별도 정규화를 넣지는 않는다.

R6의 `topm` 전체 풀 → bit/value별 stable sort → 부족한 행의 전체 fallback을 유지한다.
`argpartition`의 풀 순서는 인덱스 순서가 아니다.
그래서 풀 내부 stable sort만으로 전역적인 동점 우선순위가 정해지지는 않는다.
이번 **NEW 동점 정책**은 작은 codeword 인덱스 우선이다.
풀의 인덱스를 먼저 오름차순 정렬한다.
풀 경계의 동점이 eta번째 선택에 영향을 주면 그 행만 전체 후보로 fallback한다.

정당화는 다음과 같다. 풀 밖 점수는 풀의 최저 점수 이하이다.
선택된 마지막 점수가 그 최저 점수보다 크면 풀 밖 후보는 순서를 바꿀 수 없다.
같고 동점 후보가 풀 밖에도 있으면 전체 stable sort가 필요하다.
유효 후보가 부족한 경우에도 전체 fallback이 필요하다.
이 두 경우를 `bad`로 묶어 배치 처리한다. R6의 `for rr` 루프는 만들지 않는다.

이 정책은 R6의 기존 **동점 선택**을 바꿀 수 있다. 같은 리스트 바이트를 보장하지 않는다.
비동점 거리 순서는 보존한다. float32에서 갈라진 수학적 동점을 tolerance로 합치지 않는다.
`nfall`은 고유 중심 수가 아니라 `(중심,bit,value)` fallback 행 수다.
이번에는 경계 동점 fallback도 포함한다. loop는 codebook chunk·bit·bit-value 축에만 있다.

### 3. A1 MATLAB 단계의 수식 — [VERIFY V2]

여기서 M=2, `d=T-2`, `nd=4d`, `Q=2^B0`, `B=nd*B0`다.
T≥4를 받는다. T=3은 H의 양의 차수 조건 `T-3≥1`을 만족하지 않는다.
이 절의 W는 MATLAB의 국소 사상 변수다. 프로젝트 채널 잡음 W를 재정의하지 않는다.
상단 square root는 A_top으로 써서 eigenframe A와 구분한다.

#### 3.1 E53–57, E87–102: 격자와 Gaussian 역변환

MATLAB symbol과 Python digit의 관계는 `symbol_j=digit_j+1`이다.
확률 좌표 `x_j=lattice[digit_j]`를 선택한다.
배열 순서는 `(r의 실수,허수 교대 d개; s의 실수,허수 교대 d개)`다.
표준 실수 정규 CDF를 F_N이라 쓰면

$$
\begin{aligned}
r_\ell&=\frac{F_N^{-1}(x_{2\ell-1})+iF_N^{-1}(x_{2\ell})}{\sqrt2},\\
s_\ell&=\frac{F_N^{-1}(x_{2d+2\ell-1})+iF_N^{-1}(x_{2d+2\ell})}{\sqrt2},\quad 1\le\ell\le d.
\end{aligned}
\tag{9}
$$

E96–99의 `norminv(...,0,1/sqrt(2))`를 그대로 표현했다.
독립 **연속 균등** 확률 좌표에서만 실수·허수 분산 1/2의 Gaussian이 나온다.
유한 격자 좌표에 Gaussian 분포나 정확한 분산 1/2를 주장하지 않는다.
E129의 균등분포 설명도 유한 codebook의 분포를 뜻하지 않는다.

#### 3.2 H34–35: thnt 방정식과 벡터화 가능한 근

MATLAB H의 차수 n을 여기서는 q라 쓴다. 배치 크기 n과 구별하기 위해서다.
`theta=thnt(q,t)`가 만족하는 식은

$$
(q+1)\theta^{2q}-q\theta^{2q+2}
=1-e^{-t^2}\sum_{j=0}^{q-1}\frac{t^{2j}}{j!}=: \mathsf P_q(t^2).
\tag{10}
$$

우변은 정규화 lower incomplete gamma다.
좌변에 `z=theta²`를 대입한다.
정규화 incomplete beta의 적분을 직접 계산하면

$$
I_z(q,2)=q(q+1)\int_0^z u^{q-1}(1-u)\,du
=(q+1)z^q-qz^{q+1}.
\tag{11}
$$

미분은 `q(q+1)z^(q-1)(1-z)>0`이고 양 끝 값은 0,1이다.
따라서 내부 근은 유일하며 순방향과 역방향은

$$
\theta=\sqrt{I^{-1}_{\mathsf P_q(t^2)}(q,2)},\qquad
t=\sqrt{\mathsf P_q^{-1}\bigl(I_{\theta^2}(q,2)\bigr)}.
\tag{12}
$$

이는 H의 방정식에서 유도한 신규 수치 표현이다. 새로운 논문 수식으로 인용하지 않는다.
NumPy 배치 연산과 기존 의존성 SciPy의 beta/gamma 특수함수를 쓴다.
scalar fzero를 Python 블록 루프로 감싸지 않는다.
큰 확률은 complement를 사용해 `1-p`의 상쇄를 피한다.
MATLAB의 fzero 반복 경로나 마지막 부동소수 비트까지 복제한다는 주장은 하지 않는다.

영점은 점근식을 비교해서 처리한다.

$$
\mathsf P_q(t^2)\sim\frac{t^{2q}}{q!},\quad
I_{\theta^2}(q,2)\sim(q+1)\theta^{2q}
\quad\Longrightarrow\quad
\lim_{t\to0}\frac{\theta}{t}=[(q+1)!]^{-1/(2q)}.
\tag{13}
$$

역방향 비율은 그 역수다. `np.divide(...,where=...)`로 0/0의 선행 평가를 막는다.

#### 3.3 E107–118: 직교 분해와 radial 사상

$$
u=\frac{r^Hs}{\|r\|^2}r,\qquad \bar s=s-u,\qquad r^H\bar s=0.
\tag{14}
$$

E112와 E113–118은 각각

$$
p=\frac{r}{\|r\|}\operatorname{thnt}(T-2,\|r\|),\qquad
\bar q=\frac{\bar s}{\|\bar s\|}\operatorname{thnt}(T-3,\|\bar s\|).
\tag{15}
$$

원본도 `bar s=0`이면 `bar q=0`으로 둔다.
반면 r=0이면 E108의 투영 방향이 정의되지 않는다. 포트는 임의 방향을 만들지 않고 오류를 낸다.
짝수 Q의 대칭 등간격 격자는 1/2을 포함하지 않는다.
사용자 지정 격자가 같은 조건을 만족한다고 가정하지는 않는다.

#### 3.4 E121–134: 평행 성분과 semiunitary 행렬

$$
\begin{aligned}
v&=\begin{cases}
0,&u=0,\\
\dfrac{u}{\|u\|}\sqrt{1-e^{-\|u\|^2}}
\sqrt{(1-\|p\|^2)(1-\|\bar q\|^2)},&u\ne0,
\end{cases}\\
W&=[p,\bar q+v],\qquad
A_{\rm top}=(I_2-W^HW)^{1/2},\qquad
X=\begin{bmatrix}A_{\rm top}\\W\end{bmatrix}.
\end{aligned}
\tag{16}
$$

v는 p와 평행하고 bar q는 p에 직교한다.
a=||p||, b=||bar q||, c=||v||라 쓰면

$$
\det(I_2-W^HW)=(1-a^2)(1-b^2)-c^2
=(1-a^2)(1-b^2)e^{-\|u\|^2}>0.
\tag{17}
$$

유한 내부 입력에서는 a,b<1이다.
첫 주대각 minor도 `1-a²>0`이므로 `I₂-WᴴW`는 positive definite다.
principal Hermitian square root를 취하면

$$
X^HX=A_{\rm top}^HA_{\rm top}+W^HW
=A_{\rm top}^2+W^HW=I_2.
\tag{18}
$$

따라서 E134의 sqrtm은 Hermitian eig 재구성으로 계산한다.
Cholesky, 원소별 sqrt, 임의의 열별 정규화로 바꾸지 않는다.

#### 3.5 D54–58, D104–108: subspace와 polar gauge

Y의 상위 두 left singular vector로 `C=[C1;C2]`를 만든다.
원본의 두 번째 SVD와 gauge 제거는

$$
C_1=U_1S_1V_1^H,\qquad Q=U_1V_1^H,\qquad W=C_2Q^H.
\tag{19}
$$

무잡음, rho>0, rank(H)=2이면 `C=XR`인 unitary R이 존재한다.
이때 `C1=A_top R`, `C2=W R`이다.
A_top이 positive definite이므로

$$
\operatorname{polar}(A_{\rm top}R)
=A_{\rm top}R(R^HA_{\rm top}^2R)^{-1/2}
=A_{\rm top}RR^HA_{\rm top}^{-1}R=R.
\tag{20}
$$

따라서 `C2 Qᴴ=W`가 복원된다.
부호·위상뿐 아니라 선택된 두 차원 안의 모든 unitary basis mixing이 소거된다.
G 경로는 `YYᴴ`의 상위 두 고유벡터로 같은 subspace를 얻는다.
이 G1 계산이 전체 posterior를 subspace만으로 조건화한다는 주장은 하지 않는다.

#### 3.6 D113–130: inverse radial 단계

$$
p=W_{:,1},\quad t=W_{:,2},\quad
v=\frac{p^Ht}{\|p\|^2}p,\quad \bar q=t-v,\quad
D_0=\sqrt{(1-\|p\|^2)(1-\|\bar q\|^2)}.
\tag{21}
$$

D119–121의 근은 `a0=||v||/D0`에 대해

$$
\|u\|=\sqrt{-\log(1-a_0^2)},\qquad
u=\frac{\|u\|}{\|v\|}v.
\tag{22}
$$

D124–127은 식 (12)를 적용한다.

$$
\begin{aligned}
r&=\frac{p}{\|p\|}\sqrt{\mathsf P_{T-2}^{-1}\bigl(I_{\|p\|^2}(T-2,2)\bigr)},\\
\bar s&=\frac{\bar q}{\|\bar q\|}\sqrt{\mathsf P_{T-3}^{-1}\bigl(I_{\|\bar q\|^2}(T-3,2)\bigr)},
\qquad s=\bar s+u.
\end{aligned}
\tag{23}
$$

**원본 D121·127에는 영점 가드가 없다.**
v=0이면 u=0, bar q=0이면 bar s=0이라는 연속극한을 **NEW 확장**으로 넣었다.
v→0에서 `||u||/||v||→1/D0`이고 다른 두 비율은 식 (13)의 역수다.
원본이 이미 이 확장을 구현했다고 적지 않는다.
p=0, singular C1, D0=0에는 임의 인덱스를 만들지 않는다. 오류로 처리한다.

#### 3.7 D135–163: 정규 CDF와 hard decision

복원한 r,s의 실수 성분들을 식 (9)와 같은 순서로 z_j에 넣는다.

$$
\widehat x_j=F_N(\sqrt2z_j),\qquad
\widehat d_j=\arg\min_{0\le q<Q}(\widehat x_j-\operatorname{lattice}[q])^2.
\tag{24}
$$

D143–157의 min은 동점에서 첫 인덱스를 고른다.
NumPy argmin의 첫 인덱스 규칙은 0-기반으로 대응한다.
`np.rint`의 half-to-even이나 R6의 `floor(Q*x)`를 가져오지 않는다.
후자는 일반 alpha 격자의 최근접 규칙이 아니다.

#### 3.8 BG25–28, GB25–32, F113–116·128–130: Gray 라벨

각 좌표의 자연 이진 정수 q에 대해 Gray 정수는

$$
g=q\mathbin{\mathrm{XOR}}(q\!\gg\!1),\qquad
b_\ell=g_0\mathbin{\mathrm{XOR}}\cdots\mathbin{\mathrm{XOR}}g_\ell.
\tag{25}
$$

첫 식은 Bin2Gray, 둘째는 Gray2Bin의 prefix XOR다.
각 좌표 안에서는 MSB부터 쓴다.
F114–115의 방향은 `payload → Gray2Bin → bit2int → +1 → encoder`다.
F129–130에서는 `symbol-1 → int2bit → Bin2Gray`로 payload를 복원한다.

전역 codeword 인덱스는 R6와 같은 자리 규약으로 새로 정의한다.

$$
k=\sum_{j=0}^{nd-1}d_jQ^j,\qquad
\mathrm{bits}[k]=\mathrm{GrayBits}_{B0}(d_0)\Vert\cdots\Vert\mathrm{GrayBits}_{B0}(d_{nd-1}).
\tag{26}
$$

첫 좌표는 k의 최하위 자리지만 그 좌표의 비트 그룹은 MSB 우선이다.
`bits[k]`를 k의 단순 이진 표현으로 대신하면 안 된다.
MATLAB 인코더가 이 전역 인덱스를 정의한 것은 아니다. 식 (26)은 **NEW R6 통합 규약**이다.

### 4. MATLAB과 NumPy의 규약 차이 — [VERIFY V3]

| 지점 | 포트에서 유지하거나 명시한 것 |
|---|---|
| 1-기반 인덱스 | MATLAB symbol=digit+1. Python 내부와 최종 codeword index는 0-기반 |
| BG26·GB26의 reshape | `(B0,nd)`의 **열 하나가 좌표 하나**다. 원본의 rows 주석이 아니라 실행문을 따른다 |
| 열 우선 reshape | 한 블록의 `reshape((B0,nd),order='F').T`는 `reshape((nd,B0))`와 같다. 배치는 C-order `(n,nd,B0)`로 표현한다 |
| batch와 `(:)` | 원본 `(:)`는 블록 내부 벡터화다. 배치 전체를 Fortran reshape해 서로 다른 블록을 섞지 않는다 |
| 전치 | MATLAB `'`는 켤레전치, `.'`는 비켤레전치. 3차원 NumPy `.T`는 배치 축까지 뒤집는다. 행렬 켤레전치는 `.conj().swapaxes(-1,-2)`다 |
| SVD 반환 | MATLAB `[U,S,V]`, NumPy `(U,s,Vh)`. Q는 `U @ Vh`이며 Vh를 다시 켤레전치하지 않는다 |
| SVD 순서·크기 | 양쪽 모두 특이값 내림차순을 사용한다. D57의 full SVD는 N≥2에서 상위 두 reduced 열로 대체 가능하다 |
| N<2 | full U가 추가로 주는 열은 식별된 rank-two 신호 subspace가 아니다. G1_Y에서 거부한다 |
| eigh 순서 | NumPy는 오름차순이다. G의 마지막 두 열을 역순으로 취한다. sqrtm 재구성은 고유값과 대응 벡터를 함께 쓴다 |
| 부호·위상·중복값 | U 원소의 직접 일치를 요구하지 않는다. 선택 subspace 안의 unitary 변화는 식 (20)으로 제거된다 |
| 경계 동률 | 두 번째와 세 번째 고유값이 같으면 선택 subspace가 비유일하다. 이 경우 라이브러리 간 결정 일치를 보장하지 않는다 |
| sqrtm | Hermitian principal square root를 쓴다. Cholesky나 elementwise sqrt로 대체하지 않는다 |
| norm/reduction | 좌표 축만 줄인다. 축 없는 norm으로 batch 전체를 합치지 않는다 |
| fzero | 식 (12)·(22)의 동일 방정식 해를 구한다. 원본 반복 경로와 floating-point 비트는 보존 대상이 아니다 |
| 난수 | MATLAB과 NumPy의 같은 seed는 같은 H,Z를 뜻하지 않는다. 참조 비교에는 실제 Y를 공유한다 |

수치 guard는 double precision에서 기계 epsilon 상수배 안의 작은 이탈만 보정한다.
역함수의 반올림된 상단 1은 `nextafter(1,0)`으로 제한한다.
큰 domain 위반·비유한 값·원본 radial bracket [0,100] 밖의 해는 오류다.
이 clipping과 영점 확장은 `NEW`이며 원본 경계 동작과 같다고 주장하지 않는다.
G 경로는 Gram 형성 때문에 작은 singular value의 정보를 더 쉽게 잃는다.
따라서 Y/G의 수치 rank guard가 모든 ill-conditioned 입력에서 같은 판정을 준다는 보장은 없다.

### 5. 코드 블록 1 — 추출 대상 `code/m2_gate.py`

아래 모듈은 NumPy로 배치를 처리한다. 특수함수는 R6에도 있는 SciPy 의존성을 사용한다.
라이선스 고지까지 포함해 추출한다. 이 문서만 커밋하며 원본 R6는 수정하지 않는다.

```python
# NEW: Preserve the upstream license below when extracting this A1 port.
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

## 검산 방법

### 1. 사양 §7과의 대응

| 항목 | 실행할 검사와 범위 |
|---|---|
| 1: M=1 환원 | 동일 Cube-split codebook에서 Y·RNG 상태·κ·Phi·지표·ML 결정·오류 개수·exact LLR 비교. 지표 상대오차 <1e-6 |
| 2: Grassmann codebook | B0=1,2의 **모든** k에서 `||XᴴX-I₂||F<1e-10` |
| 3: 라벨 | 1000 무작위 인덱스 왕복과 독립 bit-group XOR 비교 |
| 4: 무잡음 복원 | G1만 검사. 2000 블록에서 Y 경로와 G 경로 모두 인덱스 완전 일치. G2는 다음 chat |
| 5: G2≤ML | G2는 범위 밖이며 이 문서로 인증하지 않음 |
| 6: 참조 구현 | Octave/MATLAB 미실행. 아래 참조 벡터 절차가 필요 |
| 7: 정밀도 | 독립 complex128 trace와 기존 float32 fast ML argmax 일치율 ≥0.999 |
| 8: 리스트 경계 | 계산된 float32 점수의 eta·eta+1 동점 횟수와 분모 기록. 전체 stable sort와 리스트를 대조 |
| 9: LDPC | 재작성하지 않음. 원본 sanity와 실제 gate 설정에서 R6 구성·rate를 확인해야 함 |

**item 1 실행:** 송신 k 생성 RNG와 채널 RNG를 분리한다.
동일 seed의 독립 RNG 세 개로 기존 channel, 새 legacy 입력, 새 `(n,T,1)` 입력을 비교한다.
출력뿐 아니라 RNG 잔여 상태를 검사하므로 추가 난수 소모도 검출한다.
같은 k에 대한 오류 **개수의 equality**로 같은 SER를 확인한다. 성능 수치를 새로 보고하지 않는다.
기존 gvec·fast ML·exact LLR을 그대로 호출한다.

**item 4 실행:** 각 B0에서 k 2000개를 먼저 고정한다.
H를 한 번 생성해 rank(H)=2를 확인하고 `Y0=sqrt(T/2)*(X@H)`를 만든다.
이는 rho=1, Z=0이다. 잡음을 더하는 `channel()`을 호출하면 무잡음 검사가 아니다.
두 G1 API의 반환값이 k와 완전히 같은지 assert한다.
하나라도 틀리거나 예외가 나면 실패다. 해당 블록 삭제·H 재추출로 분모를 바꾸지 않는다.
전체 codeword의 고정 identity-H 복원과 임의 right-unitary gauge 검사도 추가했다.
identity-H 검사는 대수적 fixture이지 폐기된 물리적 N=M 특수 경우 주장이 아니다.

추가로 일반 M의 무작위 semiunitary 입력에서 trace와 chordal score를 독립 식으로 확인한다.
작은 topm을 써서 fallback 경로를 강제로 검사한다.
H34 방정식 잔차도 확인한다. 이 잔차 검사는 MATLAB 실행 비교를 대신하지 않는다.

### 2. 코드 블록 2 — 추출 대상 `code/sanity_m2.py`

검산 코드는 모두 `NEW`다. 원본 알고리즘을 다시 구현하는 대신 앞 모듈과 R6를 호출한다.
root residual과 직접 trace만 독립 산술 대조식으로 사용한다.

```python
import json  # NEW: JSON-compatible YAML needs no additional dependency
import sys  # NEW
from pathlib import Path  # NEW
import numpy as np  # NEW
from scipy.special import gammainc  # NEW: equation-residual oracle, not a decoder
import cubesplit_gap as r6  # NEW: actual unchanged project source
import m2_gate as m2  # NEW: first Python block of GPT01


def check_m1(cfg):  # NEW: EXP section 7 item 1
    rng = np.random.default_rng(cfg['seed'])  # NEW
    cs = r6.CubeSplit(4, 1)  # NEW: small R6 fixture, not a performance sweep
    C, bits = cs.codebook()  # NEW
    k = rng.integers(cs.K, size=2000)  # NEW: fixed test inputs
    X, phi0 = C[:, k].T, r6.codebook_features(C)  # NEW
    for N in cfg['N_m1']:  # NEW: parameter loop, not a received-block loop
        for rho in cfg['rho_linear']:  # NEW: parameter loop
            a, b, c = (np.random.default_rng(cfg['seed']) for _ in range(3))  # NEW: separate identical RNGs
            Y0 = r6.channel(a, X, N, rho)  # NEW
            Y1 = m2.channel(b, X, N, rho)  # NEW
            Y2 = m2.channel(c, X[:, :, None], N, rho)  # NEW
            assert np.array_equal(Y0, Y1) and np.array_equal(Y0, Y2)  # NEW: stronger than SER equality
            assert a.bit_generator.state == b.bit_generator.state == c.bit_generator.state  # NEW
            assert m2.kappa_of(rho, cs.T) == r6.kappa_of(rho, cs.T)  # NEW
            assert m2.kappa_of(rho, cs.T, 1) == r6.kappa_of(rho, cs.T)  # NEW
            assert np.array_equal(m2.codebook_features(C), phi0)  # NEW
            phi1 = m2.codebook_features(C[:, None, :])  # NEW
            assert np.array_equal(phi1, phi0)  # NEW
            G0, G2 = r6.gram(Y0), r6.gram(Y2)  # NEW: unchanged gram
            met0, met1 = r6.ml_metrics_fast(G0, phi0), r6.ml_metrics_fast(G2, phi1)  # NEW
            rel = np.linalg.norm((met0 - met1).astype(float)) / max(np.linalg.norm(met0.astype(float)), np.finfo(float).tiny)  # NEW
            assert rel < 1e-6  # NEW: item 1
            pred0, pred1 = met0.argmax(1), met1.argmax(1)  # NEW
            assert np.array_equal(pred0, pred1)  # NEW
            assert np.count_nonzero(pred0 != k) == np.count_nonzero(pred1 != k)  # NEW: same SER, no rates reported
            L0 = r6.exact_llr_chunked(G0, phi0, bits, r6.kappa_of(rho, cs.T))  # NEW
            L1 = r6.exact_llr_chunked(G2, phi1, bits, m2.kappa_of(rho, cs.T, 1))  # NEW
            assert np.array_equal(L0, L1)  # NEW: additional unchanged-LLR integration test


def check_lists(gl, C, bits, cfg):  # NEW: EXP item 8, full stable-sort oracle
    em = max(cfg['etas']) + 1  # NEW: eta+1 needed for boundary ties
    lists, _ = m2.precompute_lists_sorted(gl, C, bits, em, cfg['chunk'], cfg['topm'])  # NEW: R6 argument order
    C32 = np.ascontiguousarray(C, dtype=np.complex64)  # NEW
    ties = {eta: 0 for eta in cfg['etas']}  # NEW: computed-score ties, not approximate ties
    for s in range(0, gl.K, cfg['chunk']):  # NEW: codebook-chunk loop
        e = min(s + cfg['chunk'], gl.K)  # NEW
        corr = m2._chordal_rows(C32, s, e)  # NEW
        for j in range(gl.B):  # NEW: bit-coordinate loop
            for b in (0, 1):  # NEW
                oracle = np.argsort(-np.where(bits[:, j] == b, corr, -1.0), axis=1, kind='stable')[:, :em]  # NEW
                assert np.array_equal(lists[s:e, j, b], oracle)  # NEW: candidate order and fallback
                selected = np.take_along_axis(corr, oracle, 1)  # NEW
                for eta in cfg['etas']:  # NEW: list-size loop
                    ties[eta] += int(np.count_nonzero(selected[:, eta - 1] == selected[:, eta]))  # NEW
    for eta in cfg['etas']:  # NEW
        print({'B': gl.B, 'eta': eta, 'boundary_equal_f32': ties[eta], 'denominator': gl.K * gl.B * 2})  # NEW: output measured counts only on execution


def check_a1(cfg):  # NEW: EXP items 2,3,4(G1),7
    rng = np.random.default_rng(cfg['seed'])  # NEW
    for B0 in (1, 2):  # NEW: T=4, B=8 and 16
        Q, alpha = 1 << B0, cfg['alpha_fixture']  # NEW: not a gate alpha choice
        lattice = alpha + np.arange(Q) * (1 - 2 * alpha) / (Q - 1)  # F91,99-100
        gl = m2.GrassLattice(4, B0, lattice)  # NEW
        C, bits = gl.codebook()  # NEW
        Xall = C.transpose(2, 0, 1)  # NEW
        err = np.linalg.norm(Xall.conj().swapaxes(1, 2) @ Xall - np.eye(2), axis=(1, 2))  # NEW
        assert np.all(err < 1e-10)  # NEW: item 2, every codeword
        k = rng.integers(gl.K, size=1000)  # NEW: item 3
        d = gl.from_index(k)  # NEW
        assert np.array_equal(gl.index(gl.from_bits(gl.bits(d))), k)  # NEW
        assert np.array_equal(gl.from_bits(bits[k]), d)  # NEW
        binary = ((d[:, :, None] >> (B0 - 1 - np.arange(B0))) & 1).astype(np.uint8)  # NEW: independent bit-group construction
        gray = np.concatenate((binary[:, :, :1], binary[:, :, 1:] ^ binary[:, :, :-1]), 2)  # BG27; NEW: independent oracle
        assert np.array_equal(gray.reshape(1000, gl.B), bits[k])  # NEW: round-trip alone misses a shared permutation
        k0 = rng.integers(gl.K, size=2000)  # NEW: item 4, fixed denominator
        X0 = gl.symbols(gl.from_index(k0))  # NEW
        for N in cfg['N_a1']:  # NEW: parameter loop
            assert N >= 2  # NEW: no rank-two recovery claim when N<2
            H = (rng.standard_normal((2000, 2, N)) + 1j * rng.standard_normal((2000, 2, N))) / np.sqrt(2)  # NEW
            assert np.all(np.linalg.matrix_rank(H) == 2)  # NEW: no redraw or deletion of failures
            Y0 = np.sqrt(gl.T / 2) * (X0 @ H)  # NEW: rho=1, Z=0; never call the noisy channel here
            assert np.array_equal(m2.decode_G1_Y(gl, Y0), k0)  # NEW: item 4 Y path
            assert np.array_equal(m2.decode_G1(gl, r6.gram(Y0)), k0)  # NEW: item 4 G path
        assert np.array_equal(m2.decode_G1_Y(gl, Xall), np.arange(gl.K))  # NEW: all-codeword identity-H fixture, not a physical special case
        R0 = rng.standard_normal((2000, 2, 2)) + 1j * rng.standard_normal((2000, 2, 2))  # NEW
        U, _, vh = np.linalg.svd(R0)  # NEW: arbitrary right-unitary mixing, not just signs
        assert np.array_equal(m2._decode_basis(gl, X0 @ (U @ vh)), k0)  # NEW: polar-gauge invariance
        ncheck = cfg['precision_blocks']  # NEW: fixed before execution
        k1 = rng.integers(gl.K, size=ncheck)  # NEW
        Y = m2.channel(rng, gl.symbols(gl.from_index(k1)), cfg['N_a1'][-1], cfg['precision_rho_linear'])  # NEW
        G, phi = r6.gram(Y), m2.codebook_features(C)  # NEW
        P = np.einsum('tmk,smk->tsk', C, C.conj())  # NEW: actual XX^H in complex128, not rounded Phi
        pred32, pred64 = np.empty(ncheck, np.int64), np.empty(ncheck, np.int64)  # NEW
        for s in range(0, ncheck, cfg['chunk']):  # NEW: batch-chunk loop
            e = min(s + cfg['chunk'], ncheck)  # NEW
            m64 = np.einsum('nts,stk->nk', G[s:e], P).real  # NEW: independent double-precision trace
            m32 = r6.ml_metrics_fast(G[s:e], phi)  # NEW
            assert np.linalg.norm(m64 - m32) / np.linalg.norm(m64) < 1e-6  # NEW
            pred32[s:e], pred64[s:e] = m32.argmax(1), m64.argmax(1)  # NEW
        assert np.mean(pred32 == pred64) >= 0.999  # NEW: item 7
        if cfg['check_lists'] and (B0 == 1 or cfg['full_lists']):  # NEW: large full setup is opt-in
            check_lists(gl, C, bits, cfg)  # NEW
        else:  # NEW: never report skipped coverage as complete
            print({'B': gl.B, 'lists': 'not run'})  # NEW


def check_generic(cfg):  # NEW: general-M feature, distance, and root-equation arithmetic
    rng = np.random.default_rng(cfg['seed'])  # NEW
    for T in (1, 4):  # NEW: dimensions of arithmetic fixtures
        for M in range(1, T + 1):  # NEW: rank-parameter loop
            A = rng.standard_normal((16, T, M)) + 1j * rng.standard_normal((16, T, M))  # NEW
            X, _ = np.linalg.qr(A, mode='reduced')  # NEW: semiunitary frames
            C = X.transpose(1, 2, 0)  # NEW
            G = r6.gram(m2.channel(rng, X, 2, 1.0))  # NEW
            direct = np.einsum('tmk,nts,smk->nk', C.conj(), G, C).real  # NEW
            fast = r6.ml_metrics_fast(G, m2.codebook_features(C))  # NEW
            assert np.linalg.norm(fast - direct) / np.linalg.norm(direct) < 1e-6  # NEW
            bits = ((np.arange(16)[:, None] >> np.arange(3, -1, -1)) & 1).astype(np.uint8)  # NEW
            C32 = np.ascontiguousarray(C, dtype=np.complex64)  # NEW
            corr = m2._chordal_rows(C32, 0, 16)  # NEW
            pair = np.einsum('tma,tlb->amlb', C.conj(), C)  # NEW: independent all-pair inner products
            assert np.allclose(corr, (np.abs(pair)**2).sum((1, 2)), rtol=1e-6, atol=1e-6)  # NEW
            lists, _ = m2.precompute_lists_sorted(None, C, bits, 4, chunk=4, topm=8)  # NEW: exercise truncated pool/fallback
            for j in range(4):  # NEW: bit-coordinate loop
                for b in (0, 1):  # NEW
                    oracle = np.argsort(-np.where(bits[:, j] == b, corr, -1.0), axis=1, kind='stable')[:, :4]  # NEW
                    assert np.array_equal(lists[:, j, b], oracle)  # NEW
    t = np.asarray([0.0, 0.001, 0.1, 0.5, 1.0, 2.0, 4.0])  # NEW: root-equation fixture, not a channel sweep
    for order in (1, 2, 3, 4):  # NEW: function-order loop
        theta = m2.thnt(order, t)  # NEW
        residual = (order + 1)*theta**(2*order) - order*theta**(2*order + 2) - gammainc(order, t*t)  # H34; NEW: residual check
        assert np.all(np.abs(residual) < 1e-12)  # NEW: arithmetic criterion, not MATLAB equivalence


if __name__ == '__main__':  # NEW
    cfg = json.loads(Path(sys.argv[1]).read_text(encoding='utf-8'))  # NEW: JSON subset of YAML
    check_m1(cfg)  # NEW
    check_generic(cfg)  # NEW
    check_a1(cfg)  # NEW
    print('completed named assertions; G2, Octave, and LDPC are not certified here')  # NEW
```

### 3. 설정과 실행

다음은 **검산 입력**이지 gate 설정 또는 측정 결과가 아니다.
`alpha_fixture=0.01`은 F58의 입력 범위에서 택했다. 최적 alpha라는 뜻이 아니다.
JSON 표기는 YAML의 부분집합이다. `configs/GPT01_sanity.yaml`로 저장하고 표준 라이브러리로 읽는다.

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

저장소 루트에서 코드 두 블록을 각각 위 경로로 추출한 후 실행한다.
아래 셸 줄은 모두 `NEW` 절차다.

```bash
PY="$HOME/miniforge3/envs/torch/bin/python"  # NEW: repository interpreter convention
mkdir -p logs  # NEW: local logs, not committed here
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$PY" code/sanity_m2.py configs/GPT01_sanity.yaml > logs/GPT01_sanity.log 2>&1  # NEW: preserve process exit status
(cd code && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$PY" sanity.py) > logs/GPT01_R6_sanity.log 2>&1  # NEW: unchanged R6/LDPC checks
```

기본 설정은 B=16 전체 리스트 setup을 실행하지 않는다. 로그에도 미실행을 표시한다.
이를 검사할 때는 `full_lists=true`로 바꾸어 별도 실행한다.
실제 gate 격자가 결정되면 그 입력으로 모든 해당 검사를 다시 수행해야 한다.
검산용 alpha를 결과를 본 뒤 gate 설정으로 승격하지 않는다.

### 4. 참조 대조 절차 — EXP §8

grassbox의 실제 Git commit을 고정하고 첨부 E/D/H/BG/GB와 비교한다.
아래 SHA-256을 upstream Git commit이라고 쓰지 않는다.
고정한 digits 200개와 실제 Y 200개를 양쪽에 동일하게 공급한다.
같은 seed의 MATLAB/NumPy 난수 일치를 가정하지 않는다.

참조 인코더에는 `digits+1`을 넣는다.
참조 X의 `(T,2,n)`과 Python `(n,T,2)`는 **축만** 옮겨 대조한다.
참조 복호기가 반환한 좌표에서는 1을 빼고 식 (26)으로 전역 인덱스를 계산한다.
codeword 오차 기준은 사양의 1e-10이며 복호 인덱스는 완전히 같아야 한다.
`.npz`에는 digits, lattice, Y, X_ref, k_ref, upstream commit을 보존한다.
원본의 영점 NaN이나 실패 사례를 삭제하지 않는다. NEW 연속극한 확장과의 차이로 기록한다.

### 5. 이 대화에서 실제 확인한 범위

첨부된 `cubesplit_gap.py` **전체를 변경 없이 import**하여 위 검산 코드를 실행했다.
M=1 환원, 일반 M trace·거리, A1 라벨·semiunitarity·무잡음 복원,
polar gauge, f32/f64 지표, 작은 codebook 리스트의 assertion이 통과했다.
이는 해당 fixture의 내부 검산이다. 독립 검산이나 실제 gate 격자 인증은 아니다.
생성된 수치 로그를 이 문서에서 새 실험 결과로 인용하지 않는다.
B=16 전체 리스트, Octave/MATLAB, G2, 서버의 원본 LDPC sanity는 실행하지 않았다.

## [VERIFY] 목록

| ID | 아직 확정하지 않은 것 | 확인 방법 |
|---|---|---|
| V1 | 공통 M 일반화와 R6 통합 | 서버에서 블록 2의 M1·일반 M 검사를 실제 R6 전체 모듈로 재실행한다. 식 (3)·(7)·(8)을 독립 대조한다 |
| V2 | A1 신규 유도와 수치적 재표현 | 식 (10)–(24), root 잔차, zero-limit를 검토한다. 고정 MATLAB 참조 벡터와 비교한다 |
| V3 | MATLAB 포트의 indexing·Gray·SVD 대응 | 공유 digits/Y로 codeword와 복호 인덱스를 비교한다. 열 우선 그룹 순서와 arbitrary unitary gauge를 별도 검사한다 |
| V4 | 실제 gate lattice/alpha | 저자가 설정 파일로 고정한 격자를 넣고 EXP §7의 2·3·4·7을 다시 실행한다 |
| V5 | 리스트 동점 정책과 B=16 전체 setup | NEW 작은-index 우선 정책을 승인한 뒤 `full_lists=true`로 전체 stable-sort oracle 및 eta 경계를 대조한다 |
| V6 | singular/near-boundary 수치 정책 | p=0·v=0·qbar=0·singular C1·작은 eigengap fixture를 구분한다. 원본 실패와 NEW 확장을 숨기지 않고 기록한다 |
| V7 | LDPC 연결과 전체 gate readiness | 변경 없는 R6 sanity, 실제 gate LDPC 설정의 구성·rate를 확인한다. A2·G2는 별도 구현·검산한다 |

[VERIFY]는 내부 fixture 통과만으로 닫지 않는다.
특히 V2·V3의 참조 비교를 하지 않은 상태에서 “원논문 복호기와 같다”고 쓰지 않는다.

## 이 결과가 바꾸는 것

`00_STATE.md` “잔여 [VERIFY]” 5번의 **“코드 없음”**에만 다음 교체안을 제안한다.
입력 사본에서는 69번째 줄이다.

> 공통 M 층과 A1 G1 코드 초안은 GPT01에 있다. 독립·참조 검산 전이다.
> A2, G2, soft 연결 및 gate 측정은 미완료다.

같은 파일 “다음 작업” 5번은 완료 처리하지 않는다.
“확정”, 기존 gate 판정, 네트워크·학습 제한을 바꾸지 않는다.
`01_PLAN.md`의 STALE 절은 인용하거나 교체하지 않는다.
`02_RESULTS_LOG.md`에 성능 결과를 추가하지 않는다. **성능은 측정 필요**다.
이 제안은 문서 안의 변경안이며 STATE/PLAN 파일을 직접 수정하지 않았다.

## ⚠️ 한계

**참조 구현과 대조하지 못했다.** Octave가 설치되어 있지 않았고 MATLAB도 실행하지 않았다.
같은 수학적 내부 사상에서 출발했지만 fzero 대신 특수함수 해를 썼다.
영점 연속극한, roundoff guard, numerical-rank 거부와 동점 선택은 NEW 변경이다.
이 차이를 숨긴 채 “그대로 재현했다”고 표현해서는 안 된다.

격자 값이 사양에서 정해지지 않아 완성된 gate codebook은 아직 하나로 고정되지 않았다.
부동소수 경계에서는 encoder contraction이 반올림되어 singular해질 수 있다.
큰 domain 위반은 clip으로 덮지 않고 오류로 처리한다.
G/Y 경로의 ill-conditioned rank 판정, 선택 subspace 경계 동률, singular polar gauge에는 동일 출력을 보장하지 않는다.

A1의 전체 K 열거는 작은 gate codebook을 위한 것이다.
int64 인덱스 범위가 허용된다는 사실은 큰 K의 메모리나 실행 가능성을 보장하지 않는다.
일반 M 코드도 모든 extreme 차원·SNR에서 수치 정확도를 인증한 것은 아니다.
SER/BLER 격차, 실행 시간, 비용 측정, 새로운 문헌은 이 산출물에 없다.

### 첨부 원본 고정 정보

아래는 첨부 바이트의 SHA-256이다. upstream Git commit의 증거가 아니다.

| 파일 | SHA-256 |
|---|---|
| GrassLatticeEncoding.m | `71567a7e569b4d692fe29d629a53c36f257a0ac5ccc1d502b0b570262e46124c` |
| GrassLatticeDecoding.m | `8df1e8b27cb3f088f7904c246f2a792f54f60bcc6652cc1c5b994a1954241d59` |
| thnt.m | `8ef5534944b14cd24a768e68cac86d3857509d080a48f47b41254d51a8badd10` |
| Bin2Gray.m | `51f7f718805f55201e9dba97c6f7b9d67fff093c4c0b7935d431a7db580dda13` |
| Gray2Bin.m | `8aa1088898d008a443cdf1d24861001af3de92abf748afe282567e55f86ab249` |
| FindAlphaOpt.m | `44cc9d49813b59bbcbcc5ee3c6af6f4096b2d1dd6e96088066d3f6acdf823688` |
| grassbox_LICENSE.txt | `6a725b2d47dec98c34a98a76239f2277acedd59e29daa45650d5c353de562ba4` |

첨부 BSD-3 고지를 첫 코드 블록에 보존했다. 코드 추출 시에도 유지한다.
