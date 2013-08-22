function [ y ] = DiffbarAbs( x, k )
%DIFFBARABS differentiable function approximating abs()
% k specifies the sharpness at x=0

y = (log(2) + log(cosh(k*x)) )/k;

end

