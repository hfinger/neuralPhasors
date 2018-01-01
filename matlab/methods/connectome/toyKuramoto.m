clear all;

N=10; %number of neurons per population
numCluster=5; % number of populations

k_intraCluster=1; %connections strenght within cluster
delay_intraClust=0.1; %connections delays within cluster

k_interCluster=1; %connections strenght between cluster
delay_interClust=0.3; %connections delays between cluster

k=1; %global scaling of connection strength;
v=1; %global conduction delay

weightedNetwork=true;
binaryNetwork=false;

f=60; %frequency of oscillators
dt=0.001; %simulation step size (in seconds)
sampling=0.001;  %saving frequency
t_max=1; % simulation time (in seconds)
sig_n=0.1;

%%
C=k_intraCluster*(ones(N,N) - diag(ones(1,N)));
D=delay_intraClust*10*(ones(N,N) - diag(ones(1,N)));

if numCluster>1
  C2 = k_interCluster*ones(numCluster*N,numCluster*N);
  D2 = delay_interClust*10*ones(numCluster*N,numCluster*N);
  for clust=1:numCluster
    ids=(clust-1)*N+1:clust*N;
    C2(ids,ids)=C;
    D2(ids,ids)=D;
  end
  C=C2;
  D=D2;
end

if weightedNetwork
%   C=C.*randn(size(C));
  C=C.*gamrnd(5,1,size(C));
end
if binaryNetwork
  C=(C > rand(size(C)));
end


%%
startPhase = 2*pi*rand(size(C,1),1, 'double');
startState = bsxfun(@plus, startPhase, (0:dt:1)*2*pi*f);

Y = runKuramoto(C,D,k,f,v,t_max,dt,sampling,sig_n,d,verbose,approx,invertSin,startState );

%%
globPhase = angle(sum(exp(1i*Y),1));
localPhase = mod(bsxfun(@minus,Y,globPhase),2*pi);
clusterPhase = zeros(numCluster,size(localPhase,2));
clusterOrderCoeff = zeros(numCluster,size(localPhase,2));
for clust=1:numCluster
  ids=(clust-1)*N+1:clust*N;
  clusterPhase(clust,:) = angle(sum(exp(1i*Y(ids,:)),1)); 
  clusterOrderCoeff(clust,:) = abs(sum(exp(1i*Y(ids,:)),1)/size(Y,1)); 
end


%%
figure(1); clf;
imagesc(sin(Y))

figure(2); clf;
imagesc(1000*diff(Y,[],2)/(2*pi));
xlabel('time [ms]')

colormap jet;
colorbar;

% figure(3); clf;
% imagesc(localPhase); 
% colormap hsv; 
% colorbar

figure(4); clf;
imagesc(C);
xlabel('source neuron')
ylabel('target neuron')

% figure(5); clf;
% imagesc(clusterPhase)

% figure(6); clf;
% plotPolarPhases( localPhase(:,1), numCluster, N );

figure(8);
plot(clusterOrderCoeff')
xlabel('time [ms]')
clusterNames=num2cell(num2str((1:size(clusterOrderCoeff,1))'));
legend(clusterNames{:})
