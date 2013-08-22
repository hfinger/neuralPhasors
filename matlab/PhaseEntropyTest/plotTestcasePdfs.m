clear all;
close all;
clc;
addpath('include')

%% A) vary num clusters, std in clusters=0, background is uniform, vary clusterNeurons/bgNeurons
stdClusters=0;
numClusters=[1 4 7];
numNeurons=[100 500 900];
densityEstSigma=0.05;
for i=1:length(numNeurons)
    for j=1:length(numClusters)
        genVarsBg = genClusters( 1000-numNeurons(i), 1 , 0 ); %equidistant background
        genVars = [genClusters( numClusters(j), round(numNeurons(i)/numClusters(j)) , stdClusters ); genVarsBg];
        phase = genPhase( genVars );
        evalOptFcn(phase,densityEstSigma,['../../Documentation/images/testcaseA/NumNeurons' num2str(numNeurons(i)) 'NumClusters' num2str(numClusters(j))]);
    end
end


%% B) vary num clusters, vary std in clusters, no background
stdClustersLn=[-3 -1 1];
stdClusters=exp(stdClustersLn);
numClusters=[1 4 7];
numNeurons=1000;
densityEstSigma=0.05;
genVarsBg = []; %no additional background
for i=1:length(stdClusters)
    for j=1:length(numClusters)
        genVars = [genClusters( numClusters(j), round(numNeurons/numClusters(j)) , stdClusters(i) ); genVarsBg];
        phase = genPhase( genVars );
        evalOptFcn(phase,densityEstSigma,['../../Documentation/images/testcaseB/stdCluster' num2str(stdClustersLn(i)) 'NumClusters' num2str(numClusters(j))]);
    end
end

