function [ genVars ] = genClusters( numClusters, numNeurons, stdDev )
%GENCLUSTERS Summary of this function goes here
%   Detailed explanation goes here

if numNeurons>0 && numClusters>0
    genVars=[numNeurons*ones(numClusters,1) ...
        linspace(0,2*pi-2*pi/numClusters,numClusters)' ...
        stdDev*ones(numClusters,1)];
else
    genVars=zeros(0,3);
end

