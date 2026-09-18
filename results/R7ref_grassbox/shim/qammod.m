function y = qammod(x, M, varargin)
% qammod — MATLAB Communications Toolbox 호출 규약 대체 (정사각 QAM 한정).
%
% Octave의 communications 패키지 qammod 는 'UnitAveragePower' 옵션을 받지 않는다.
% grassbox의 ExpMapEncoding.m 은 그 옵션을 쓰므로, 같은 Gray 성상도에 평균전력
% 정규화만 더한 이 wrapper를 path 앞에 둔다. grassbox 파일은 수정하지 않는다.
%
% 성상도 자체가 forge qammod 와 일치하는지는 gen_ref.m 이 실행 전에 대조한다.

  uap = false;
  i = 1;
  while i <= numel(varargin)
    if ischar(varargin{i}) && strcmpi(varargin{i}, 'UnitAveragePower')
      uap = logical(varargin{i+1}); i = i + 2;
    else
      i = i + 1;
    end
  end

  k = log2(M);
  if mod(k,2) ~= 0, error('qammod shim: square QAM only (M = 4^n)'); end
  h = k/2; L = 2^h;

  xi = double(x(:));
  gi = bitshift(xi, -h);        % 상위 h비트 -> I축 Gray 인덱스
  gq = bitand(xi, L-1);         % 하위 h비트 -> Q축 Gray 인덱스
  bi = gray2bin_idx(gi, h);
  bq = gray2bin_idx(gq, h);

  I = 2*bi - (L-1);             % -(L-1) .. (L-1), 왼쪽에서 오른쪽
  Qv = (L-1) - 2*bq;            % +(L-1) .. -(L-1), 위에서 아래
  y = I + 1i*Qv;

  if uap, y = y / sqrt(2*(M-1)/3); end
  y = reshape(y, size(x));
end

function b = gray2bin_idx(g, n)
  b = g;
  sh = 1;
  while sh < n
    b = bitxor(b, bitshift(b, -sh));
    sh = sh * 2;
  end
end
