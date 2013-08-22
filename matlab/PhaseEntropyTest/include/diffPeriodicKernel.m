function [ K ] = diffPeriodicKernel( phi, sigma )
%DIFFPERIODICKERNEL derivative of "von Mises" distribution

kappa=1/sigma^2;
const=2*pi*besseli(0, kappa);
K = -sin(phi).*exp(kappa*cos(phi)) / const;


end

