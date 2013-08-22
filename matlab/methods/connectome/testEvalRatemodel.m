
load('datasimu.mat');

C = bsxfun(@rdivide,SC,sum(SC,2));
% C=SC;
D=D;
k=0.8;
tau=20;
v=10;
t_max=500;
dt=0.0001;
sampling=10;
sig_n=0.25;

d=0;
verbose=true;
t_rm = 20;

% CONNECTION MATRIX =================================
R = length(C);

% START SIMULATION ================================
disp('beginning dynamics ...');
tic;
Y = runRatemodel(C,D,k,tau,v,t_max,dt,sampling,sig_n,d,verbose);
save('rawSimSignalK0p8','Y');
if ~isempty(find(~isfinite(Y),1)), 
    error('simulation failed !!!'); 
end
time_dyna = toc;

% COMPUTE BOLD SIGNAL ==============================
% using the nonlinear balloon-windkessel model...
disp('beginning bold calculation ...');
tic;
Ybold = mBOLDs(Y,dt*sampling,t_rm);
time_bold = toc;

% COMPUTE SOME BASIC BOLD SIGNAL ANALYSES ==========
% settings for bold averages
xsec = 2000; xgap = 500;        % in msec, window size and spacing

T = size(Ybold,2)*dt*sampling*1000;
t0 = 1:xgap:T-xsec+1; t0 = ceil(t0/(dt*sampling*1000));
te = xsec:xgap:T; te = ceil(te/(dt*sampling*1000));

% compute bold averages
Ybold_w = zeros(R,length(t0));
for w=1:length(t0), Ybold_w(:,w) = mean(Ybold(:,t0(w):te(w)),2); end

% remove NaNs, get average bold signal over whole brain, and regress out
Ybold_w(isnan(Ybold_w)) = 0;
Ybold_w_mean = mean(Ybold_w,1);
Ybold_w_reg = zeros(R,length(t0));
if ~isempty(Ybold_w), for i=1:R, [~,~,Ybold_w_reg(i,:)] = regress(Ybold_w(i,:)',Ybold_w_mean'); end, end

if nargout==0, 
    save('testEvalRatemodelResults','Ybold_w_reg','time_dyna','time_bold','Ybold_w');
end

disp('... all done ...');
