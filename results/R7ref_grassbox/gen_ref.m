% gen_ref.m — R7 참조 벡터 생성 (grassbox, M=2, T=4)
%
% EXP_M2_gate_spec.md §8. Python 포트 대조용 참조 벡터를 만든다.
% grassbox의 .m 파일은 읽기만 한다. 수정하지 않는다.
%
% 실행: conda run -n octave octave-cli gen_ref.m
%
% 대상
%   A1 Grass-Lattice, M=2, T=4, B'in{1,2}  -> B = 8B' bit
%   A2 Exp-Map,       M=2, T=4, Q in{4,16} -> B = 4*log2(Q) bit
%
% 규약은 Testscript_Structured_M2.m 을 그대로 따른다.
%   H     = (randn(M,N)+1i*randn(M,N))/sqrt(2)
%   nv    = M/(T*snr),  Noise = sqrt(nv/2)*(randn(T,N)+1i*randn(T,N))
%   Y     = X*H + Noise          (송신 X에 별도 스케일 없음. 잡음 분산이 SNR을 담는다)
% grassbox에는 채널 생성 함수가 없으므로 위 코드를 여기에 그대로 둔다.

1;  % script file

% ---------- MATLAB Communications Toolbox 대체 (라벨 규약용) ----------
% int2bit/bit2int 는 툴박스 함수라 Octave에 없다. MSB-first 규약 그대로 구현한다.
function b = i2b(x, n)   % 정수 -> n비트 열벡터 (MSB first)
  b = zeros(n,1);
  for i = 1:n
    b(i) = bitget(x, n-i+1);
  end
end

function x = b2i(b, n)   % 비트열 -> n비트씩 묶어 정수 (MSB first)
  L = numel(b)/n;
  Bm = reshape(b(:), n, L);
  w  = 2.^(n-1:-1:0);
  x  = (w * Bm)';
end

% ---------- 설정 ----------
GRASSBOX = '/home/HTJ/ng/code/ref/grassbox/StructuredConstellations/functions';
SHIM     = fullfile(pwd, 'shim');
addpath(GRASSBOX);
pkg load statistics communications;

% MATLAB의 qammod/qamdemod 'UnitAveragePower' 옵션이 Octave communications 1.2.7 에
% 없다. shim/ 의 wrapper를 path 앞에 두되, 성상도가 forge 구현과 같은지 먼저 대조한다.
forge_map = {};
for qq = [4 16]
  forge_map{end+1} = qammod((0:qq-1)', qq);
end
addpath(SHIM, '-begin');
jj = 0;
for qq = [4 16]
  jj = jj + 1;
  d = max(abs(qammod((0:qq-1)', qq) - forge_map{jj}));
  printf('shim qammod vs forge, Q=%2d : maxdiff = %g\n', qq, d);
  if d ~= 0, error('qammod shim does not match the communications package'); end
end

T = 4; M = 2; N = 4;
NREF   = 200;          % 인덱스 0..199
SNRdB  = 15;           % 참조 Y 생성 SNR
SEED   = 20260919;
snr    = 10^(SNRdB/10);
NoiseVar = M/(T*snr);

% Grass-Lattice alpha: Testscript_Structured_M2.m 은 T=4,M=2,B'=2 에 0.15 를 쓴다
% ([3, Table 1]). B'=1 의 표값은 공개 코드에 없다. 두 경우 모두 0.15 로 고정하고
% 그 사실을 여기에 기록한다. 참조 벡터의 목적은 포트 대조이므로 alpha 값 자체는
% Python 쪽과 같기만 하면 된다.
ALPHA_GL = 0.15;

printf('=== R7 ref vectors (grassbox) ===\n');
printf('T=%d M=%d N=%d  NREF=%d  SNRdB=%g  NoiseVar=%.17g  seed=%d\n', ...
       T, M, N, NREF, SNRdB, NoiseVar, SEED);
printf('alpha (Grass-Lattice) = %.17g\n', ALPHA_GL);
printf('octave %s\n', version());
fflush(stdout);

% ---------- A1 Grass-Lattice ----------
for Bp = [1 2]
  randn('state', SEED); rand('state', SEED);   % 설정마다 같은 자리에서 시작

  P       = 2^Bp;
  lattice = ALPHA_GL + (0:P-1)*(1-2*ALPHA_GL)./(P-1);
  Nsym    = 2*M*(T-M);          % 좌표 수 = 8
  Nbits   = Nsym*Bp;            % B = 8B'
  K       = 2^Nbits;

  printf('\n--- A1 Grass-Lattice  B''=%d  P=%d  Nbits=%d  K=%d ---\n', Bp, P, Nbits, K);
  printf('lattice = '); printf('%.17g ', lattice); printf('\n');

  Xall = zeros(T, M, NREF);
  Yall = zeros(T, N, NREF);
  Hall = zeros(M, N, NREF);
  tx_sym = zeros(Nsym, NREF);
  rx_sym = zeros(Nsym, NREF);
  tx_idx = (0:NREF-1)';
  rx_idx = zeros(NREF,1);
  nfail_noiseless = 0;

  for ii = 1:NREF
    k = tx_idx(ii);

    % 인덱스 -> 비트 -> (Gray2Bin) -> 좌표 심볼   [Testscript_Structured_M2.m]
    tx_bits    = i2b(k, Nbits);
    tx_bitsbin = Gray2Bin(tx_bits, Bp);
    s          = b2i(tx_bitsbin, Bp) + 1;        % 1..P, 길이 Nsym
    X          = GrassLatticeEncoding(M, lattice, s);

    H     = (randn(M,N) + 1i*randn(M,N))/sqrt(2);
    Noise = sqrt(NoiseVar/2)*randn(T,N) + 1i*sqrt(NoiseVar/2)*randn(T,N);
    Y     = X*H + Noise;

    r          = GrassLatticeDecoding(M, lattice, Y);
    rx_bitsbin = reshape(cell2mat(arrayfun(@(v) i2b(v-1,Bp), r(:)', 'UniformOutput', false)), [], 1);
    rx_bits    = Bin2Gray(rx_bitsbin, Bp);
    kh         = 0;
    for j = 1:Nbits, kh = kh*2 + rx_bits(j); end

    % 무잡음 왕복 확인 (복호기·인코더 일치 여부)
    r0 = GrassLatticeDecoding(M, lattice, X*H);
    if any(r0(:) ~= s(:)), nfail_noiseless = nfail_noiseless + 1; end

    Xall(:,:,ii) = X; Yall(:,:,ii) = Y; Hall(:,:,ii) = H;
    tx_sym(:,ii) = s; rx_sym(:,ii) = r; rx_idx(ii) = kh;
  end

  nerr = sum(rx_idx ~= tx_idx);
  printf('noiseless roundtrip failures : %d / %d\n', nfail_noiseless, NREF);
  printf('index errors at %g dB        : %d / %d\n', SNRdB, nerr, NREF);

  label_convention = ['A1 Grass-Lattice: index k -> %d bits MSB-first -> Gray2Bin(.,B'') ' ...
                      '-> B''-bit groups MSB-first -> coordinate symbol 1..P. ' ...
                      'Decoder returns symbols 1..P; invert with int2bit(s-1,B'') then Bin2Gray.'];
  label_convention = sprintf(label_convention, Nbits);

  fname = sprintf('ref_A1_T4_Bp%d.mat', Bp);
  save('-v7', fname, 'T','M','N','Bp','P','lattice','Nsym','Nbits','K', ...
       'SNRdB','NoiseVar','SEED','ALPHA_GL', ...
       'Xall','Yall','Hall','tx_idx','rx_idx','tx_sym','rx_sym', ...
       'nerr','nfail_noiseless','label_convention');
  printf('saved %s\n', fname);
  fflush(stdout);
end

% ---------- A2 Exp-Map ----------
for Q = [4 16]
  randn('state', SEED); rand('state', SEED);

  Bq    = log2(Q);              % bits/QAM symbol
  Nsym  = M*(T-M);              % QAM 심볼 수 = 4
  Nbits = Nsym*Bq;
  K     = 2^Nbits;

  printf('\n--- A2 Exp-Map  Q=%d  Nbits=%d  K=%d ---\n', Q, Nbits, K);

  Xall = zeros(T, M, NREF);
  Yall = zeros(T, N, NREF);
  Hall = zeros(M, N, NREF);
  tx_sym = zeros(Nsym, NREF);
  rx_sym = zeros(Nsym, NREF);
  tx_idx = (0:NREF-1)';
  rx_idx = zeros(NREF,1);
  nfail_noiseless = 0;

  for ii = 1:NREF
    k = tx_idx(ii);

    tx_bits = i2b(k, Nbits);
    s       = b2i(tx_bits, Bq);           % 0..Q-1, Gray 변환 없음 (qammod 내부 Gray)
    X       = ExpMapEncoding(M, Q, s);

    H     = (randn(M,N) + 1i*randn(M,N))/sqrt(2);
    Noise = sqrt(NoiseVar/2)*randn(T,N) + 1i*sqrt(NoiseVar/2)*randn(T,N);
    Y     = X*H + Noise;

    r       = ExpMapDecoding(M, Q, Y);
    rx_bits = reshape(cell2mat(arrayfun(@(v) i2b(v,Bq), r(:)', 'UniformOutput', false)), [], 1);
    kh = 0;
    for j = 1:Nbits, kh = kh*2 + rx_bits(j); end

    r0 = ExpMapDecoding(M, Q, X*H);
    if any(r0(:) ~= s(:)), nfail_noiseless = nfail_noiseless + 1; end

    Xall(:,:,ii) = X; Yall(:,:,ii) = Y; Hall(:,:,ii) = H;
    tx_sym(:,ii) = s; rx_sym(:,ii) = r; rx_idx(ii) = kh;
  end

  nerr = sum(rx_idx ~= tx_idx);
  printf('noiseless roundtrip failures : %d / %d\n', nfail_noiseless, NREF);
  printf('index errors at %g dB        : %d / %d\n', SNRdB, nerr, NREF);

  % QAM 성상도 자체를 기록으로 남긴다 (라벨 규약의 실체)
  qam_map = qammod((0:Q-1)', Q, 'UnitAveragePower', true);
  printf('qam_map (symbol 0..%d):\n', Q-1);
  for q = 1:Q, printf('  %2d  % .17g % +.17gi\n', q-1, real(qam_map(q)), imag(qam_map(q))); end

  if T == 4, alpha_em = 0.3; else, alpha_em = NaN; end
  label_convention = sprintf(['A2 Exp-Map: index k -> %d bits MSB-first -> %d-bit groups MSB-first ' ...
                              '-> QAM symbol 0..%d, no external Gray step (qammod Gray mapping, ' ...
                              'UnitAveragePower). qam_map is stored. alpha=%g (T=4).'], ...
                             Nbits, Bq, Q-1, alpha_em);

  fname = sprintf('ref_A2_T4_Q%d.mat', Q);
  save('-v7', fname, 'T','M','N','Q','Bq','Nsym','Nbits','K', ...
       'SNRdB','NoiseVar','SEED','alpha_em','qam_map', ...
       'Xall','Yall','Hall','tx_idx','rx_idx','tx_sym','rx_sym', ...
       'nerr','nfail_noiseless','label_convention');
  printf('saved %s\n', fname);
  fflush(stdout);
end

printf('\n=== done ===\n');
