
function [cost, grad, cc] = froNorm(SC, FC, gSC, lambda, gamma, k)

b = inv(eye(size(gSC))-k*gSC);
bb = b*b';

% calculate cost-function
cstr = norm(gSC(:)-SC(:))^2;
cfnc = norm(FC(:) -bb(:))^2;
cc.c1 = cstr;
cc.c2 = cfnc;
cost = lambda*cc.c1 + gamma*cc.c2;

% calculate gradient matrix
if nargout > 1
  reg = 2*(gSC-SC);
  
  phi = FC - bb;                                                            % phi == outer derivative
  tmp1 = b' * phi * bb';
  tmp2 = b' * phi' * bb;
  des = 2*(-k)*(tmp1 + tmp2);
  
  grad = lambda*reg + gamma*des;
  grad(logical(eye(size(grad)))) = 0;                                       % remove diagonal entries
else
  grad = zeros(size(gSC));
end

end