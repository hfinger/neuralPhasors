function [out] = evalTestcaseSurf(testcase)
addpath('include')

if strcmp(testcase,'A')
    %% A) vary num clusters, std in clusters=0, background is uniform, vary clusterNeurons/bgNeurons
    stdClusters=0;
    numClusters=1:9;
    numNeuronsBg=0:100:1000;
    densityEstSigma=0.05;
    for i=1:length(numNeuronsBg)
        disp(num2str(i))
        for j=1:length(numClusters)
            genVarsUniform = genClusters( numNeuronsBg(i), 1 , 0 ); %equidistant background
            genVarsCluster = genClusters( numClusters(j), round((1000-numNeuronsBg(i))/numClusters(j)) , stdClusters );
            genVars = [genVarsCluster; genVarsUniform];
            phase = genPhase( genVars );
            outvars(i,j) = evalOptFcn(phase,densityEstSigma);
        end
    end
    
    out.xData = numNeuronsBg/10;
    out.xLabel = 'uniform distributed [\%]';
    out.yData = numClusters;
    out.yLabel = 'number of clusters';
elseif strcmp(testcase,'B')
    %% B) vary num clusters, vary std in clusters, no background
    stdClustersLn=-5:0.5:1;
    stdClusters=exp(stdClustersLn);
    numClusters=1:9;
    numNeuronsBg=1000;
    densityEstSigma=0.05;
    genVarsUniform = []; %no additional background
    for i=1:length(stdClusters)
        disp(num2str(i))
        for j=1:length(numClusters)
            genVars = [genClusters( numClusters(j), round(numNeuronsBg/numClusters(j)) , stdClusters(i) ); genVarsUniform];
            phase = genPhase( genVars );
            outvars(i,j) = evalOptFcn(phase,densityEstSigma);
        end
    end
    
    out.xData = stdClustersLn;
    out.xLabel = 'ln(std in clusters)';
    out.yData = numClusters;
    out.yLabel = 'number of clusters';
elseif strcmp(testcase,'C')
    %% C) vary position of cluster, vary fraction in background
    stdClusters=0;
    numNeuronsBg=1000:-100:100;
    posSecondCluster=(9*pi/10):-(pi/10):0;
    numClusters=2;
    densityEstSigma=0.05;
    for i=1:length(numNeuronsBg)
        disp(num2str(i))
        for j=1:length(posSecondCluster)
            genVarsUniform = genClusters( numNeuronsBg(i), 1 , 0 ); %equidistant background
            genVarsCluster=[round((1000-numNeuronsBg(i))/numClusters), 0, stdClusters;...
                round((1000-numNeuronsBg(i))/numClusters), posSecondCluster(j), stdClusters];
            
            genVars = [genVarsCluster; genVarsUniform];
            phase = genPhase( genVars );
            outvars(i,j) = evalOptFcn(phase,densityEstSigma);
        end
    end
    
    out.xData = numNeuronsBg/10;
    out.xLabel = 'uniform distributed [\%]';
    out.yData = posSecondCluster;
    out.yLabel = 'cluster-distance';
end

fnames=fieldnames(outvars);
for i=1:length(fnames)
    out.data.(fnames{i}) = reshape([outvars.(fnames{i})],size(outvars,1),size(outvars,2));
end