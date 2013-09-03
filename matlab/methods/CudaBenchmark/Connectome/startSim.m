numRois = 660;
t_max=1;
approx = false;
useGPU = true;

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
phase = runKuramoto(C,D,t_max,approx,useGPU);
disp(toc)

