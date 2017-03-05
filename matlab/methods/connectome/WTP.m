function [ m ] = WTP( u, e0, u0, r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

m = 2*e0./(1 + exp(r*(u0 - u)));

end

