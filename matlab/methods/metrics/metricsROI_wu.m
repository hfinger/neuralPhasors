function   [perROI,lable] = metricsROI_wu(SC)

paths = dataPaths();
load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));

ut = triu(true(size(SC)),+1)+triu(true(size(SC)),+1)';                      % de-select main diagonal
inter = logical(zeros(size(SC)));                                           % mask for interhemispheric connections
inter(1:length(SC)/2,length(SC)/2+1:end) = logical(1);
inter = logical(inter + inter');
intra = logical(logical(ut - inter) + logical(ut - inter)');                % mask for intrahemispheric connections

%% measures per ROI: 66x1

perROI.degree = sum((SC~=0))';                                              % total degrees (undirected)
perROI.intradegree = sum((SC.*intra) ~= 0,2);                               % degree per ROI of only intrahemispheric connections
perROI.interdegree = sum((SC.*inter) ~= 0,2);                               % degree per ROI of only interhemispheric connections
perROI.strengths = sum(SC+SC',2);                                           % total weights, symm: in- & out (see van den Heuvel & Sporns 2011, Fig. 2a)

[in,os,~] = strengths_dir(SC);                                              % since undirected: in = os = strengths_und(SC)
perROI.outStrength = os';
perROI.inpStrength = in';                                                   % post-RowNormalization: inputs cdf to 1 (graph no longer symmetric!)

perROI.intraStrength = sum(SC.*intra,2);
perROI.interStrength = sum(SC.*inter,2);


perROI.clustCoef = clustering_coef_wu(SC);                                  % see Fagiolo (2007) Phys Rev E 76:026107
perROI.betwCentr = betweenness_wei(1./SC);                                  % input should be a distance matrix, give 1/SC
perROI.roiSize = avg_roi_size;

perROI.shortestPaths = sum(distance_wei(1./SC),1)';                         % get all-pairs shortest paths, calc. average per ROI


%% measures of hemispheric correlation: 33x1

hLeft = SC(1:length(SC)/2,1:length(SC)/2);
hRight= SC(length(SC)/2+1:length(SC),length(SC)/2+1:length(SC));

% cosine similarity [0,1] for SC>0
perROI.cosSimiC = (sum(hLeft.*hRight,1)./(sqrt(sum(hLeft.^2,1)).*sqrt(sum(hRight.^2,1))))'; % cosSimi along columns (outputs)
perROI.cosSimiR = sum(hLeft.*hRight,2)./(sqrt(sum(hLeft.^2,2)).*sqrt(sum(hRight.^2,2)));    % cosSimi along rows (inputs)                                                                  %

perROI.jacSimi = (sum(min(hLeft,hRight),1) ./ sum(max(hLeft,hRight),1))';   % generalized Jaccard similarity
perROI.corSimi = corr(hLeft,hRight);                                        % Pearson's linear correlation coeff
perROI.corSimi = perROI.corSimi(eye(length(SC)/2)==1);                      % ALL ROI highly +correlated, useless?

perROI.MIdSimi = matching_ind_und(SC>0)';                                   % matching index (on binary graph), req. sparseness
% MId is symmetric, hence MIdSimi(h1) = MIdSimi(h2)

%% communities and modularity
[perROI.Ci, perROI.Q] = modularity_und(SC,.75);                             % see Newman (2006)
perROI.pIdx = participation_coef(SC,perROI.Ci,0);                           % sum(pIdx(1:33)) < sum(pIdx(33:66))

% testwise calculation of communities per hemisphere
[aaa1,bbb1] = modularity_und(SC(1:length(SC)/2,1:length(SC)/2),.25);
aaa1=vertcat(aaa1,zeros(33,1));


lable = [cellstr('degree');cellstr('strengths');cellstr('outStrength');cellstr('inpStrength');
         cellstr('intraStrength');cellstr('interStrength');cellstr('clustCoef');cellstr('betwCentr');cellstr('roiSize');...
         cellstr('cosSimiC');cellstr('cosSimiR');cellstr('jacSimi');cellstr('corSimi');cellstr('MIdSimi');...
         cellstr('Ci');cellstr('Q');cellstr('pIdx'); cellstr('shortestPaths')];