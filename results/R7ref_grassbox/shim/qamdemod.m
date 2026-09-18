function x = qamdemod(y, M, varargin)
% qamdemod — qammod.m shim의 역. 정사각 QAM 한정, 좌표별 최근접 슬라이싱.

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
  if mod(k,2) ~= 0, error('qamdemod shim: square QAM only (M = 4^n)'); end
  h = k/2; L = 2^h;

  yv = y(:);
  if uap, yv = yv * sqrt(2*(M-1)/3); end

  bi = round((real(yv) + (L-1))/2);
  bq = round(((L-1) - imag(yv))/2);
  bi = min(max(bi, 0), L-1);
  bq = min(max(bq, 0), L-1);

  gi = bitxor(bi, bitshift(bi, -1));
  gq = bitxor(bq, bitshift(bq, -1));
  x = bitshift(gi, h) + gq;
  x = reshape(x, size(y));
end
