function [ phase ] = genPhase( genVars )
%GENPHASE generate phase variables for given clusters
%   each cluster consists of 
%   genVars=[numNeurons, meanPhase, stdPhase; ...]


%% generate phase variables:
tmpPhase=cell(1,size(genVars,2));
for j=1:size(genVars,1)
    numNeurons=genVars(j,1);
    meanPhase=genVars(j,2);
    stdPhase=genVars(j,3);
    tmpPhase{j}=mod(meanPhase+stdPhase*randn(numNeurons,1),2*pi);
end
phase=cell2mat(tmpPhase');

end

