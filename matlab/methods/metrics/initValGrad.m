function   initValGrad(SC, FC)
% function to probe the error landscape of the objective functions in the
% gradient learning procedure

%% get data

paths = dataPaths();
%load('/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIds.mat') % get vector that sorts ROIs, eventually replace this
resortIds = 1:66;

if (nargin<2)
  % load SC with sp = 0.5 , without added homotopic connections
  % load(fullfile(paths.workdir,'pebel','20150414_SAR_Metrics','ConnectomeMetrics','1SC.mat'));
  % SC = hSC;
  
  % load SC with sp = 0.0 (original data), without added homotopic connections
  load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));
  SC = avg_ci;
  SC(isnan(SC)) = 0;
  SC = SC + SC';
  SC = normGraph(SC, avg_roi_size, 'ROIprd', false, 0);
  %SC = bsxfun(@rdivide,SC,sum(SC,2));
  
  dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150125_controls_fs_funconn_bponetrial_lcmv_hilbert_3_30_entmirrored.mat']);
  tmp = mean(dataEegTmp.coh_all([1:13 15 17:20],:,5:6,:,:,7),3);
  tmp = mean(tmp,1);
  FC = squeeze(tmp);
end

%% define masks

h1 = logical(diag(ones(length(SC)/2,1),+length(SC)/2));                   % get masks for LH2RH and
h2 = logical(diag(ones(length(SC)/2,1),-length(SC)/2));                   % RH2LH homotopic connections
h = h1 + h2;

nodiag = not(eye(size(SC)));                                                % deselect main diagonal
ut = triu(true(size(SC)),+1);                                               % select upper triangular matrix
sut = (size(SC,1)*(size(SC,1)-1))./2;                                       % number of elements in upper triangular matrix
mut = (1:numel(SC))';                                                       % mask that can be set to upper triangular matrix

inter = logical(zeros(size(SC)));                                           % mask for interhemispheric connections
inter(1:length(SC)/2,length(SC)/2+1:end) = logical(1);
intra = logical(ut - inter); 

inter = logical(inter + inter');
intra = logical(intra + intra');

% outer for loop for different noise distributions
SC = bsxfun(@rdivide,SC,sum(SC,2));
noise = rand(size(SC));
noise(logical(eye(size(SC)))) = 0;
noise = bsxfun(@rdivide,noise,sum(noise,2));
for m = 0:0.05:1
    nSC = m*SC + (1-m)*noise;
    [gSC, deltaSC, dev, glob] = gradDesc(nSC, FC, 3, .65, 0, m);
    close all hidden;
    disp(strcat('m=',num2str(m),' glob. perf=', num2str(glob.CORRGC)))
end

% fFC = bsxfun(@rdivide,FC,sum(FC,2));
% corr(SC(nodiag),fFC(nodiag));
% corr(gSC(intra),fFC(intra));
% corr(gSC(inter),fFC(inter));

%%
idx = 1;
for i = 0:0.05:1
load(strcat('/work/pebel/Results/initialCondition/extra',num2str(i),'/gSC.mat'))
A(:,idx) = cell2mat({glob(:,1)}); idx = idx + 1;
end

figure();plot(A(1,:))
