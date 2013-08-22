function [ C, D ] = genNetwork( params )
%GENNETWORK Summary of this function goes here
%   Detailed explanation goes here
weightedNetwork=params.weightedNetwork;
binaryNetwork=params.binaryNetwork;
N=params.N;
k_intraCluster=params.k_intraCluster;
delay_intraClust=params.delay_intraClust;
numCluster=params.numCluster;
k_interCluster=params.k_interCluster;
delay_interClust=params.delay_interClust;



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


end

