function [ c ] = corrTrig( mat1, mat2 )
%CORRTRIG Summary of this function goes here
%   Detailed explanation goes here

ids = find(triu(ones(size(mat1)),1));

c = corr( mat1(ids), mat2(ids));

end

