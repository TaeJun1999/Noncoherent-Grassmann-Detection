1;
function b = i2b(x,n), b=zeros(n,1); for i=1:n, b(i)=bitget(x,n-i+1); end, end
function x = b2i(b,n), L=numel(b)/n; Bm=reshape(b(:),n,L); w=2.^(n-1:-1:0); x=(w*Bm)'; end
addpath('/home/HTJ/ng/code/ref/grassbox/StructuredConstellations/functions');
pkg load statistics;
T=4;M=2;N=4;Bp=1;P=2;a=0.15;lat=a+(0:P-1)*(1-2*a)/(P-1);
randn('state',20260919); rand('state',20260919);
for k=0:199
  s=b2i(Gray2Bin(i2b(k,8),Bp),Bp)+1;
  X=GrassLatticeEncoding(M,lat,s);
  H=(randn(M,N)+1i*randn(M,N))/sqrt(2);
  r=GrassLatticeDecoding(M,lat,X*H);
  if any(r(:)~=s(:))
    printf('k=%3d  tx=[%s]  rx=[%s]  wrongpos=[%s]\n', k, num2str(s'), num2str(r'), num2str(find(r(:)~=s(:))'));
  end
end
printf('diag done\n');
