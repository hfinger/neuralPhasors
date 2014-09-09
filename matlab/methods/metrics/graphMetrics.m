function   metr = graphMetrics(SC, nSamples, path)

% compare metrics to results reported in Table 1, p15778 
% van den Heuvel and Sporns (2011), 'Rich-Club Organization'

% graphMetrics should ONLY analyze and plot,
% NOT do any changes (eg normalization, sparseness) by itself

if nargin < 3
    path = pwd;
end

%% metrics per ROI: 66 x 1

metr.perROI = metricsROI_wu(SC);

%% metrics per connection (ROI1, ROI2): 66 x 66

metr.perConn = metricsConn_wu(SC);

maxOOM = max(SC);                                                           % get max entry per row
minOOM = min(SC+eye(size(SC)));                                             % get min entry per row (ignore main diagonal)

%% global network metrics: 1x1

metr.perGlob = metricsGlobal_wu(SC, nSamples);


%% complete analysis & generate some plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LH->LH  %%%  RH->LH %%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%% LH->RH  %%%  RH->RH %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

imagesc(SC)                                                               
colorbar
print(strcat(path,'/SC'),'-dpng')

imagesc(metr.perConn.euclDist)                                                               
colorbar
print(strcat(path,'/euclDist'),'-dpng')                                     % plot indicates closeness of homologous ROI

hist(metr.perROI.degree,15)
print(strcat(path,'/totalDeg'),'-dpng')

hist(metr.perROI.clustCoef,15)
print(strcat(path,'/clustCoef'),'-dpng')


%% rejected analysis:

% 1) modularization: destroying LH/RH matrix structure & homotopy access
% see Honey et al (2007) for reference

[Ci1, Q1] = modularity_dir(SC(1:33,1:33),.5);                               % what's the effect of introducing homo-
%[Ci1 Q1]  = modularity_und(SC(1:33,1:33),.5);                              % topic connections on modularity?
[On, Wr] = reorder_mod(SC(1:33,1:33),Ci1);                                  % see also reorderMAT.m

imagesc(Wr)                                                                
colorbar 
print(strcat(path,'/SCmod'),'-dpng')