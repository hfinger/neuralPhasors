function [ f ] = densityEstimate( phi0, phi, sigma )
%DENSITYESTIMATE Summary of this function goes here
%   Detailed explanation goes here

f = sum( periodicKernel( phi0-phi,sigma ) ) / length(phi);

end

