rng(taskId)

%%
paths = dataPaths();
resortIds = 1:66;

load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));
empSC = avg_ci;
empSC(isnan(empSC)) = 0;
empSC = empSC + empSC';
empSC = normGraph(empSC, avg_roi_size, 'ROIprd', true, 0);

dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150125_controls_fs_funconn_bponetrial_lcmv_hilbert_3_30_entmirrored.mat']);
tmp = mean(dataEegTmp.coh_all([1:13 15 17:20],:,5:6,:,:,7),3);
tmp = mean(tmp,1);
empFC = squeeze(tmp);
clear dataEegTmp;

triuIds = find(triu(ones(66,66),1));


%% define parameters:
p.useFrob = false;
p.useRowRenorm = 'no';
p.constrainPos = false;
p.k = 0.65;
p.lambda = 0;
p.gamma = 1;
p.negPenalty = 0.01;
p.maxSteps = 65536;
p.savepath = 'noRenorm_noConstrain_noReg_moreNoiseInterp_runLonger_unittest';
p.unittest = true;
p.useRMSprop = true;
p.rmsPropDecay = 0.9;
p.initLearnRate = 1e-4;
p.saveAtIters = [0 2.^(0:16)];
p.breakAtCost = 0.999999;

mkdir(p.savepath);

%% unittest:
% p.unittest = true;
% initSC = empSC;
% [learnSC, dev] = optimizeSC(initSC, empSC, empFC, p);

%% optimize empSC:
% p.unittest = false;
% initSC = empSC;
% [learnSC, dev, logLearnSC] = optimizeSC(initSC, empSC, empFC, p);

% [simFCcov, simFCcor] = sar(logLearnSC{1,1}{end}, 0.65);
% corr( simFCcov(triuIds), empFC(triuIds))

%% outer for loop for different noise distributions

numRepeat = 1;
noiseAmount = 0:0.05:1;
logLearnSC = cell(numRepeat, length(noiseAmount));
finalLearnSC = cell(numRepeat, length(noiseAmount));
dev = cell(numRepeat, length(noiseAmount));
for j=1:numRepeat
  disp(['j=' num2str(j)])
  noise = rand(size(empSC));
  noise(logical(eye(size(empSC)))) = 0;
  noise = bsxfun(@rdivide,noise,sum(noise,2));
  for m = 1:length(noiseAmount)
    disp(['m=' num2str(m)])
    initSC = (1-noiseAmount(m))*empSC + noiseAmount(m)*noise;
    [finalLearnSC{j,m}, dev{j,m}, logLearnSC{j,m}] = optimizeSC(initSC, empSC, empFC, p);
  end
end

%%
save(fullfile(p.savepath, ['taskId' num2str(taskId) '.mat']));
