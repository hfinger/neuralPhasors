function   Corr = postEval(gSC, k, SC, FC , SCMetr)
% ggSC = repmat(k2,1,66).*gSC;

%% load data, set parameters
paths = dataPaths();
% load('/work/pebel/Results/gradDesc/reg0 (normeq)/gSC.mat')

if (nargin<3)
  % load default SC and calculate metrics of SC
  load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));
  load(fullfile(paths.databases,'SC_Bastian','eeg_20150125_controls_fs_funconn_bponetrial_lcmv_hilbert_3_30_entmirrored.mat'));
  SC = avg_ci;
  SC(isnan(SC)) = 0;
  SC = SC + SC';
  SC = normGraph(SC, avg_roi_size, 'ROIprd', false, 0);
  
  SCMetr = graphMetrics(SC, 0);                                             % nsteps = 0, no need to evaluate global metrics
end

%% define some masks
h1 = logical(diag(ones(length(gSC)/2,1),+length(gSC)/2));                   % get masks for LH2RH and
h2 = logical(diag(ones(length(gSC)/2,1),-length(gSC)/2));                   % RH2LH homotopic connections
h = h1 + h2;

ut = triu(true(size(gSC)),+1);                                              % select upper triangular matrix
sut = (size(gSC,1)*(size(gSC,1)-1))./2;                                     % number of elements in upper triangular matrix
mut = (1:numel(gSC))';                                                      % mask that can be set to upper triangular matrix

inter = logical(zeros(size(gSC)));                                          % mask for interhemispheric connections
inter(1:length(gSC)/2,length(gSC)/2+1:end) = logical(1);
intra = logical(ut - inter);                                                % mask for intrahemispheric connections

%% calculate correlations
SC = bsxfun(@rdivide, SC, sum(SC,2));                                       % calc. directed diff. between row-norm. matrices
dSC = gSC - SC;

disp(strcat('corr(sum(SC,1)''',',sum(dSC,1)''',') = ',num2str(corr(sum(SC,1)',sum(dSC,1)'))))

%% metrics per ROI: 66 x 1
Corr.perROI.degree.all1 = corr(SCMetr.perROI.degree(:), sum(dSC,1)');
Corr.perROI.strengths.all1 = corr(SCMetr.perROI.strengths(:), sum(dSC,1)');  % !
Corr.perROI.clustCoef.all1 = corr(SCMetr.perROI.clustCoef(:), sum(dSC,1)');  % !
Corr.perROI.betwCentr.all1 = corr(SCMetr.perROI.betwCentr(:), sum(dSC,1)');  % !
Corr.perROI.roiSize.all1 = corr(SCMetr.perROI.roiSize(:), sum(dSC,1)');      % !

% removed from analysis: sum(dSC,2) is approx. 0 due to normalized SC, gSC
% Corr.perROI.degree.all2 = corr(SCMetr.perROI.degree(:), sum(dSC,2));
% Corr.perROI.strengths.all2 = corr(SCMetr.perROI.strengths(:), sum(dSC,2));
% Corr.perROI.clustCoef.all2 = corr(SCMetr.perROI.clustCoef(:), sum(dSC,2));
% Corr.perROI.betwCentr.all2 = corr(SCMetr.perROI.betwCentr(:), sum(dSC,2));
% Corr.perROI.roiSize.all2 = corr(SCMetr.perROI.roiSize(:), sum(dSC,2));

%% replicate ROI metrics along a given dimension: 66 x 66
% repmat(SCMetr.perROI.degree(:),66,1);                                     % horz = repmat(SCMetr.perROI.degree(:),1,66); horz = horz(:);
% kron(SCMetr.perROI.degree, ones(66,1));                                   % vert = repmat(SCMetr.perROI.degree(:)',66,1); vert = vert(:);

Corr.perDup.degree.all1 = corr(kron(SCMetr.perROI.degree, ones(66,1)), dSC(:));
Corr.perDup.strengths.all1 = corr(kron(SCMetr.perROI.strengths, ones(66,1)), dSC(:));
Corr.perDup.clustCoef.all1 = corr(kron(SCMetr.perROI.clustCoef, ones(66,1)), dSC(:));
Corr.perDup.betwCentr.all1 = corr(kron(SCMetr.perROI.betwCentr, ones(66,1)), dSC(:));
Corr.perDup.roiSize.all1 = corr(kron(SCMetr.perROI.roiSize, ones(66,1)), dSC(:));


% removed from analysis: correlation is approx. 0 due to normalized SC, gSC
% Corr.perDup.degree.all2 = corr(repmat(SCMetr.perROI.degree(:),66,1), dSC(:));
% Corr.perDup.strengths.all2 = corr(repmat(SCMetr.perROI.strengths,66,1), dSC(:));
% ...

%% metrics per connection: 66 x 66
Corr.perConn.euclDist.all = corr(SCMetr.perConn.euclDist(:), dSC(:));
Corr.perConn.euclDist.inter = corr(SCMetr.perConn.euclDist(inter), dSC(inter));
Corr.perConn.euclDist.intra = corr(SCMetr.perConn.euclDist(intra), dSC(intra));

Corr.perConn.floydDist.all = corr(SCMetr.perConn.floydDist(:), dSC(:));
Corr.perConn.floydDist.inter = corr(SCMetr.perConn.floydDist(inter), dSC(inter));
Corr.perConn.floydDist.intra = corr(SCMetr.perConn.floydDist(intra), dSC(intra));

Corr.perConn.floydNoEdges.all = corr(SCMetr.perConn.floydNoEdges(:), dSC(:));
Corr.perConn.floydNoEdges.inter = corr(SCMetr.perConn.floydNoEdges(inter), dSC(inter));
Corr.perConn.floydNoEdges.intra = corr(SCMetr.perConn.floydNoEdges(intra), dSC(intra));

Corr.perConn.EBC.all = corr(SCMetr.perConn.EBC(:), dSC(:));                 % !
Corr.perConn.EBC.inter = corr(SCMetr.perConn.EBC(inter), dSC(inter));
Corr.perConn.EBC.intra = corr(SCMetr.perConn.EBC(intra), dSC(intra));       % !

%% stepwise linear model

lmROI = LinearModel.stepwise(horzcat(SCMetr.perROI.degree(:), SCMetr.perROI.strengths(:), ...
                                     SCMetr.perROI.clustCoef(:),SCMetr.perROI.betwCentr(:), ...
                                     SCMetr.perROI.roiSize(:)),...
                             sum(dSC,1)', 'Upper' , 'linear');
                         
coeffs = table2array(lmROI.Coefficients(:,1)); % 5 terms: intercept 1, x1, x3, x4, x3:x4
% when comparing coefficients, remember that variable ranges are not normalized!
% so coeffs is not an indicator of weighted importance, please see:
% ranges = table2array(lms.VariableInfo(:,2))
linvar = lmROI.Formula.Terms;

lmConn = LinearModel.stepwise(horzcat(SCMetr.perConn.euclDist(:),SCMetr.perConn.floydDist(:), ...
                                      SCMetr.perConn.floydNoEdges(:),SCMetr.perConn.EBC(:), ...
                                      SCMetr.perROI.MIdSimi(:)),...    
                             dSC(:), 'Upper' , 'linear');
 
nodiag = ~eye(size(SC));                                                    % to remove diagonal entries in metrices      
nodiag = nodiag(:);                                                         % - is not effective at all

lmDegree= kron(SCMetr.perROI.degree, ones(66,1)); lmDegree = lmDegree(nodiag);
lmStrength = kron(SCMetr.perROI.strengths, ones(66,1)); lmStrength = lmStrength(nodiag);
lmClustCoef= kron(SCMetr.perROI.clustCoef, ones(66,1)); lmClustCoef= lmClustCoef(nodiag);
lmbetwCentr= kron(SCMetr.perROI.betwCentr, ones(66,1)); lmbetwCentr= lmbetwCentr(nodiag);
lmRoiSize = kron(SCMetr.perROI.roiSize, ones(66,1)); lmRoiSize= lmRoiSize(nodiag);

lmEuclDist = SCMetr.perConn.euclDist(nodiag); 
lmFloydDist = SCMetr.perConn.floydDist(nodiag); 
lmFloydEdge = SCMetr.perConn.floydNoEdges(nodiag);
lmEBC = SCMetr.perConn.EBC(nodiag);
lmMIdSimi = SCMetr.perROI.MIdSimi(nodiag);


lmMix = LinearModel.stepwise(horzcat(SC(nodiag), h(nodiag), ...
                                     lmDegree, lmStrength, lmClustCoef, lmbetwCentr, lmRoiSize, ...
                                     lmEuclDist, lmFloydDist, lmFloydEdge, lmEBC, lmMIdSimi), ...
                             dSC(nodiag), 'Upper' , 'linear');
                           
% obtain metric-based estimate of changes to SC
rdSC = zeros(size(dSC));
rdSC(nodiag) = lmMix.Fitted;                                               
rdSC(logical(eye(size(rdSC)))) = 0;
rgSC = SC + rdSC;
                           
% lmneg = LinearModel.stepwise(horzcat(SC(:), h(:), ...
%                                      nodiag.*kron(SCMetr.perROI.degree, ones(66,1)), nodiag.*kron(SCMetr.perROI.strengths, ones(66,1)), ...
%                                      nodiag.*kron(SCMetr.perROI.clustCoef, ones(66,1)), nodiag.*kron(SCMetr.perROI.betwCentr, ones(66,1)), ...
%                                      nodiag.*kron(SCMetr.perROI.roiSize, ones(66,1)), ...
%                                      nodiag.*SCMetr.perConn.euclDist(:),nodiag.*SCMetr.perConn.floydDist(:), ...
%                                      nodiag.*SCMetr.perConn.floydNoEdges(:),nodiag.*SCMetr.perConn.EBC(:), ...
%                                      nodiag.*SCMetr.perROI.MIdSimi(:)),...    
%                              reshape(dSC<0,numel(dSC),1).*dSC(:), 'Upper' , 'linear');     
%                            
% rneg = reshape(lmneg.Fitted, size(dSC));                                    % estimate of changes to SC
% rneg(logical(eye(size(rneg)))) = 0;
%                            
% lmpos = LinearModel.stepwise(horzcat(SC(:), h(:), ...
%                                      nodiag.*kron(SCMetr.perROI.degree, ones(66,1)), nodiag.*kron(SCMetr.perROI.strengths, ones(66,1)), ...
%                                      nodiag.*kron(SCMetr.perROI.clustCoef, ones(66,1)), nodiag.*kron(SCMetr.perROI.betwCentr, ones(66,1)), ...
%                                      nodiag.*kron(SCMetr.perROI.roiSize, ones(66,1)), ...
%                                      nodiag.*SCMetr.perConn.euclDist(:),nodiag.*SCMetr.perConn.floydDist(:), ...
%                                      nodiag.*SCMetr.perConn.floydNoEdges(:),nodiag.*SCMetr.perConn.EBC(:), ...
%                                      nodiag.*SCMetr.perROI.MIdSimi(:)),...    
%                              reshape(dSC>0,numel(dSC),1).*dSC(:), 'Upper' , 'linear'); 
%                            
% rpos = reshape(lmpos.Fitted, size(dSC));                                    % estimate of changes to SC
% rpos(logical(eye(size(rpos)))) = 0;
% rdSC = rpos+rneg;                                                           % corr(rdSC(:),dSC(:))                                                          
% rgSC = SC + rdSC;
                           

%rgSC = bsxfun(@rdivide,rgSC,sum(rgSC,2));

for kk = 0.4:0.05:1
  [covGC, corGC] = sar(rgSC,kk);                                             % acquire SAR(rgSC)
  
  tGC = triu(covGC,1);
  tFC = triu(FC   ,1);
  
  % (1) global prediction error (correlation)                                 % (over upper triangular matrices)
  glob.CORRGC = corr(tGC(ut), tFC(ut));                                       % this is unequal dev(end)
  
  glob.INTERGC = corr(tGC(inter), tFC(inter));                                % on inter- or intrahemispherical connections
  glob.INTRAGC = corr(tGC(intra), tFC(intra));
  
end

figure();imagesc(dSC);colorbar();colormap(b2r(min(dSC(:)),max(dSC(:))));
figure();imagesc(rdSC);colorbar(); colormap(b2r(min(dSC(:)),max(dSC(:))));

%% correlate local model errors with metrics 
 
% [covGC, corGC] = sar(gSC,k);                                                % acquire SAR(rgSC)
% 
% tGC = triu(covGC,1);
% tFC = triu(FC   ,1);
% 
% % (1) global prediction error (correlation)                                 % (over upper triangular matrices)
% glob.CORRGC = corr(tGC(ut), tFC(ut));                                       % this is unequal dev(end)
% glob.INTERGC = corr(tGC(inter), tFC(inter));                                % on inter- or intrahemispherical connections
% glob.INTRAGC = corr(tGC(intra), tFC(intra));
% 
% 
% % (2) local prediction error
% [lGC, pGC, dGC] = fit_2D_data(FC(ut), corGC(ut),'no');
% eGC = reshape(dGC, size(FC(ut)));
% 
% figure()
% up = triu(ones(length(SC)),1);
% up(~~up)= eGC;
% imagesc(up); colorbar(); colormap(b2r(min(eGC),max(eGC)));
% title('local errors of SAR(gSC)');
% 
% perROI.degree, ones(66,1)), nodiag.*kron(SCMetr.perROI.strengths, ones(66,1)), ...
%                                      nodiag.*kron(SCMetr.perROI.clustCoef, ones(66,1)), nodiag.*kron(SCMetr.perROI.betwCentr, ones(66,1)), ...
%                                      nodiag.*kron(SCMetr.perROI.roiSize, ones(66,1)), .
%                                    
% % correlations on 66 x 1 metrics & DIRECTED model error
% % - relate this to Fig. 'local errors of SAR(gSC)': ROIs are systematically over- or underestimated
% corr(SCMetr.perROI.strengths, sum(up, 1)')
% corr(SCMetr.perROI.clustCoef, sum(up, 1)')
% corr(SCMetr.perROI.betwCentr, sum(up, 1)')
% corr(SCMetr.perROI.roiSize, sum(up, 1)')
% corr(mean(SCMetr.perConn.floydDist,2), sum(up, 1)')
% 
% % correlations on 66 x 66 metrics & DIRECTED model error
% corr(SCMetr.perConn.euclDist(ut),up(ut))                                    % this decreases, compared to Holger's
% corr(SCMetr.perConn.floydDist(ut),up(ut))                                   % !
% corr(SCMetr.perConn.floydNoEdges(ut),up(ut))
% corr(SCMetr.perConn.EBC(ut),up(ut))
% corr(SCMetr.perROI.MIdSimi(ut),up(ut))
%   
% % do this also for asymmetry measures

%% systematic interpolation between SC, gSC

idx = 1;
for a = 0:0.01:1
  interSC = a.*gSC + (1-a)*SC;
  %for kk = 0.4:0.1:0.4
    [covGC, corGC] = sar(interSC,k);                                       % acquire SAR(rgSC)
    
    tGC = triu(covGC,1);
    tFC = triu(FC   ,1);
    inter(idx,1) = corr(tGC(ut), tFC(ut));                                 
    idx = idx +1;
  %end
end

% disp(sum([0;inter] < [inter;1]) - 1)                                       % strictly monotonically increasing?

%% obtain global performance (correlation) for gSC of variable lambda

nodiag = ~eye(size(SC));                                                    % to remove diagonal entries in metrices      
nodiag = nodiag(:);                                                         % - is not effective at all

lmDegree= kron(SCMetr.perROI.degree, ones(66,1)); lmDegree = lmDegree(nodiag);
lmStrength = kron(SCMetr.perROI.strengths, ones(66,1)); lmStrength = lmStrength(nodiag);
lmClustCoef= kron(SCMetr.perROI.clustCoef, ones(66,1)); lmClustCoef= lmClustCoef(nodiag);
lmbetwCentr= kron(SCMetr.perROI.betwCentr, ones(66,1)); lmbetwCentr= lmbetwCentr(nodiag);
lmRoiSize = kron(SCMetr.perROI.roiSize, ones(66,1)); lmRoiSize= lmRoiSize(nodiag);

lmEuclDist = SCMetr.perConn.euclDist(nodiag); 
lmFloydDist = SCMetr.perConn.floydDist(nodiag); 
lmFloydEdge = SCMetr.perConn.floydNoEdges(nodiag);
lmEBC = SCMetr.perConn.EBC(nodiag);
lmMIdSimi = SCMetr.perROI.MIdSimi(nodiag);

idx = 1;
for reg = 0%:0.05:1.0
  [~,dSC, ~, ~] = gradDesc(SC, FC, 3, .65, reg);
  close all;
  % IF ALREADY CALCULATED, THEN: load instead
  
  lmMix = LinearModel.stepwise(horzcat(SC(nodiag), h(nodiag), ...
                                       lmDegree, lmStrength, lmClustCoef, lmbetwCentr, lmRoiSize, ...
                                       lmEuclDist, lmFloydDist, lmFloydEdge, lmEBC, lmMIdSimi), ...
                               dSC(nodiag), 'Upper' , 'linear');
  
  % obtain metric-based estimate of changes to SC
  rdSC = zeros(size(dSC));
  rdSC(nodiag) = lmMix.Fitted;
  rdSC(logical(eye(size(rdSC)))) = 0;
  rgSC = SC + rdSC;
  
  idy = 1;
  for kk = 0.1:0.05:1
    %[covGC, ~] = sar(rgSC,kk);                                              % acquire SAR(rgSC)
    b = inv(eye(size(rgSC))-kk*rgSC);                                        %
    bb = b*b';                                                               %
    covGC = bb;                                                              %
    autocov = sqrt(diag(covGC));                                             %
    tGC = triu(covGC,1);
    tFC = triu(FC   ,1);
    
    % (1) global prediction error (correlation)                             % (over upper triangular matrices)
    CORRGC(idy,1) = corr(tGC(ut), tFC(ut));                                 % this is unequal dev(end)
    idy = idy+1;
    %     glob.INTERGC = corr(tGC(inter), tFC(inter));                      % on inter- or intrahemispherical connections
    %     glob.INTRAGC = corr(tGC(intra), tFC(intra));
  end
  CORRELATION(idx,1) = max(CORRGC);
  idx = idx + 1;
  
end



%% do plotting and saving

path = strcat(paths.localTempDir,'/Results/','postEval');
if ~exist(path,'dir')
  mkdir(path);
end

savefilename = fullfile(path,'/');
save([savefilename 'correlations.mat'],'gSC','SC', 'dSC', 'SCMetr','Corr', 'lmMix');

end

%%
function [cov, cor] = sar(gSC, k) 
if isscalar(k)
  b = inv(eye(size(gSC))-k*gSC);
else
  b = inv(eye(size(gSC))-repmat(k,1,66).*gSC);
end
bb = b*b';

cov = bb;
autocov = sqrt(diag(cov));
cor = cov ./ (autocov * autocov');
end