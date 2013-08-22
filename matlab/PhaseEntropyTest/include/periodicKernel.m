function [ K ] = periodicKernel( phi, sigma )
%PERIODICKERNEL von Mises distribution

kappa=1/sigma^2;
const=2*pi*besseli(0, kappa);
K = exp(kappa*cos(phi)) / const;


end

