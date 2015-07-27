function   perROI = metricsROI_wu(CIJ)

load('SC4Holger_avg_roi_size.mat'); 

perROI.degree = sum((CIJ~=0)); %define degree of each node                  % get total degrees (undirected)
perROI.strenths = sum(CIJ+CIJ',2);                                          % total weights, symm: in- & out (see van den Heuvel & Sporns 2011, Fig. 2a)
[in,os,~] = strengths_dir(CIJ);                                             % since undirected: in = os = strengths_und(CIJ)
perROI.outStrength = os';
perROI.inpStrength = in';                                                   % post-RowNormalization: inputs cdf to 1 (graph no longer symmetric!)

perROI.clustCoef = clustering_coef_wu(CIJ);                                 % see Fagiolo (2007) Phys Rev E 76:026107
perROI.betwCentr = betweenness_wei(1./CIJ);                                 % input should be a distance matrix, give 1/CIJ
perROI.roiSize = avg_roi_size;

%% communities and modularity
[perROI.Ci, perROI.Q] = modularity_und(CIJ,.75);                            % see Newman (2006)
perROI.pIdx = participation_coef(CIJ,perROI.Ci,0);                          % sum(pIdx(1:33)) < sum(pIdx(33:66))

% testwise calculation of communities per hemisphere
[aaa1,bbb1] = modularity_und(CIJ(1:length(CIJ)/2,1:length(CIJ)/2),.25);  
aaa1=vertcat(aaa1,zeros(33,1));