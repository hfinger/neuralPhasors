function [cov, cor] = sar(gSC, k) 
if isscalar(k)
  b = inv(eye(size(gSC))-k*gSC);
else
  b = inv(eye(size(gSC))-repmat(k,1,66).*gSC);
end
bb = b*b';

cov = bb;
autocov = sqrt(diag(cov));
cor = cov ./ (autocov * autocov');
end