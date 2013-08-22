function [outvars] = evalOptFcn( phase, densityEstSigma, plotdir )
%EVALPHASE Summary of this function goes here
%   Detailed explanation goes here

addpath('include')
n=length(phase);

%% coherence:
outvars.coherence1 = abs( sum( exp(1i*phase) ) ) / n;
outvars.coherence2 = abs( sum( exp(2i*phase) ) ) / n;
outvars.coherence3 = abs( sum( exp(3i*phase) ) ) / n;
outvars.coherence2minus1 = outvars.coherence2 - outvars.coherence1;
outvars.coherence3minus1 = outvars.coherence3 - outvars.coherence1;

%% resubstitution estimate of entropy:
densityEval=zeros(size(phase));
for j=1:n
    densityEval(j) = log2(densityEstimate( phase(j) , phase, densityEstSigma ));
end
outvars.entropy=-sum(densityEval)/n;
outvars.mirroredEntropy=DiffbarAbs(outvars.entropy+0.5,10);

%% entropyMod = integral -f(x) (1 - P(x, optEntropy, 0.2]) Log(f(x))
densityEvalMod=zeros(size(phase));
for j=1:n
    density=densityEstimate( phase(j) , phase, densityEstSigma );
    densityEvalMod(j) = log2(density) * 4200.09/(2911.28 + exp(5*density));
end
outvars.modifiedEntropy=-sum(densityEvalMod)/n;

%% circular variance:
totalPhase=sum(exp(1i*phase))/n;
outvars.asymmetry=totalPhase*conj(totalPhase);
outvars.circVar=1-sqrt( DiffbarAbs(outvars.asymmetry,10) );

%% circular sparseness:
%measured by pairwise distance compared to 4/pi (= theoretic
%pairwise distance of uniform distribution)
% outvars.circSparse=0;
% for j=1:length(phase)
%     tmp1 = exp(1i*phase(j))+exp(1i*phase);
%     outvars.circSparse=outvars.circSparse+mean(DiffbarAbs(sqrt(DiffbarAbs(tmp1.*conj(tmp1),10))-4/pi,10));
% end
% outvars.circSparse=outvars.circSparse/length(phase);

%% mean derivative:
% outvars.meanDiff=0;
% for j=1:n
%     outvars.meanDiff = outvars.meanDiff+DiffbarAbs(sum( diffPeriodicKernel( phase(j)-phase,densityEstSigma ) ),10) / n;
% end
% outvars.meanDiff=outvars.meanDiff/n;

%% mean derivative squared:
% outvars.meanDiffSquare=0;
% for j=1:n
%     outvars.meanDiffSquare = outvars.meanDiffSquare+sum( diffPeriodicKernel( phase(j)-phase,densityEstSigma ) )^2 / n;
% end
% outvars.meanDiffSquare=outvars.meanDiffSquare/n;

%% mean derivative weighted:
% outvars.meanDiffW=0;
% for j=1:n
%     outvars.meanDiffW = outvars.meanDiffW+ (DiffbarAbs(sum( diffPeriodicKernel( phase(j)-phase,densityEstSigma ) ),10) / n) /densityEstimate( phase(j) , phase, densityEstSigma );
% end
% outvars.meanDiffW=outvars.meanDiffW/n;

%% number of neurons within each cluster:
outvars.clusterSize=0;
outvars.clusterSize2=0;
for j=1:n
    %tmp1 = (exp(1i*phase(j))+exp(1i*phase))/2;
    %dist = sqrt(DiffbarAbs(tmp1.*conj(tmp1),10));
    %distAngle = 2*acos(dist);
    %outvars.clusterSize = outvars.clusterSize + mean(periodicKernel( distAngle, 0.2 ));
    outvars.clusterSize = outvars.clusterSize + mean(periodicKernel( phase(j)-phase, 0.05 )); % the same but without diffbarAbs
    outvars.clusterSize2 = outvars.clusterSize2 + mean(periodicKernel( phase(j)-phase, 0.05 )).^2; % the same but without diffbarAbs
end
outvars.clusterSize=outvars.clusterSize/n;
outvars.logClusterSize=log(outvars.clusterSize/n);
outvars.numClust=1/outvars.clusterSize;
outvars.clusterSize2=outvars.clusterSize2/n;
outvars.numClust2=1/outvars.clusterSize2;


%the same but faster computation:
% pairwiseDist=phase*ones(size(phase))'-(phase*ones(size(phase))')';
% outvars.clusterSize=mean(periodicKernel( pairwiseDist(:), 0.05 ));

%combine with entropy:
outvars.clusterSizeAndEntropy=real(1.7*outvars.entropy+outvars.clusterSize);

%% plots:
if nargin > 2
    %plotPolarPhases( phase, plotdir );
    plotPolarKernelDensity( phase, plotdir, densityEstSigma );
    plotKernelDensity( phase, plotdir, densityEstSigma );
end

end

