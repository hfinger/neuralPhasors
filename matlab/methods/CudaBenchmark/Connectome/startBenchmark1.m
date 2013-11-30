function executionTime = startBenchmark1( useGpu )
%% A simple benchmark running Kuramoto simulation...

if nargin<1
  useGpu = false;
end

t_max = 3;

numRois = 660;
approx = false;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

C=abs(randn(numRois)).^2;
D=abs(150+100*randn(numRois));

C(isnan(C)) = 0;
C = C + C';

D(isnan(D)) = 0;
D = D + D';

C = bsxfun(@rdivide,C,sum(C,2));

tic;
phase = runKuramoto(C,D,t_max,approx,useGpu);

executionTime = toc;
disp(['execution time: ' num2str(executionTime)])
