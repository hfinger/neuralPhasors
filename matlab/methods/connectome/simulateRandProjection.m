clear all;

N_I = 2; % number of sources
N_O = 5; % number of electrodes
T = 10000; % number of simulated timesteps

W = randn(N_O,N_I); % some random mixing matrix

I = randn(N_I,T); % some random source signal
I = bsxfun(@minus, I, mean(I,1)); % standardize source signal
I = bsxfun(@rdivide, I, std(I,[],1));

O = W * I; % electrode signal is given by the mix of the source signals

disp('------ Statistics ------')

corW = corr(W');
corO = corr(O');
disp(['theoretically calculated corr(W)=' num2str(corW(1,2)) ' corresponds to simulated corr(output)=' num2str(corO(1,2)) ])

covW = bsxfun(@minus,W,mean(W,2))*bsxfun(@minus,W,mean(W,2))' / size(W,2); %cov(W');
covO = bsxfun(@minus,O,mean(O,2))*bsxfun(@minus,O,mean(O,2))' / size(O,2); %cov(O');
disp(['theoretically calculated cov(W)=' num2str(covW(1,2)) ' corresponds to simulated cov(output)/N_I=' num2str(covO(1,2)/N_I) ])
