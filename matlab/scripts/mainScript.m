% mainScript to organize calling of all other functions

clear all; clc; close all;

% addScriptPaths(); % run('/net/store/nbp/phasesim/src_pebel/matlab/addScriptPaths.m')
% dataPaths();

load('dti_20141209_preprocessed.mat');                                      % get SC matrix, use average SC

% recalculate avg_ci, only take connections present in >= 75% of subjects to remove outliers?
% our matrices have way more connections than data used in van den Heuvel & Sporns (2011)
% resulting in differences of metrics. --> threshold connections/

% mat = zeros(66,66);
% for i = 1:22
%     if (isempty(ci{1,i})) continue; end
%     ci{1,i}(isnan(ci{1,i}))=0;
%     mat = mat + ci{1,i};
%  
% end
% mat = mat/(22-3);


SC = avg_ci; 
SC(isnan(SC)) = 0;
SC = SC + SC';                                                              % make SC symmetric --
                                                                            % symm.con.: ROI normalization
                                                                            % non-symm.con.: row normalization, sparsification
% trigIds = find(tril(ones([66 66]),0));


%%

sparse = 0.5:0.1:0.6;  % 1:1:8;                                             % [0.6, ...., 0.75[, graph should not decompose
heuristics = 1:1:2; % 1:1:7;
hScale = 2; % 12
kScale = 12;


SCNorm = cell(length(sparse),1);      
SCMetr = cell(length(sparse),1);   
hSC    = cell(length(heuristics),hScale,length(sparse));
hMetr  = cell(length(heuristics),hScale,length(sparse));

for sp = 1:length(sparse)                                                   % iterate over network sparseness
    
    path = strcat('Results\sp',num2str(sparse(sp)));
    mkdir(path);
    
    % tic; SPRS = metricsGlobal_wu(SC06, nSamples); toc;                      graph permutation is very expensive
    % Elapsed time is 273.425889 seconds.                                     for densely connected networks
    % tic; FULL = metricsGlobal_wu(SC, nSamples); toc;                        0.7 vs 0.0 sparseness
    % Elapsed time is 1857.124529 seconds.
    
    SCNorm{sp,1} = normGraph(SC, avg_roi_size, 'ROIsum', false, sparse(sp));
    SCMetr{sp,1} = graphMetrics(SCNorm{sp,1}, path);
    
    hSC(:,:,sp) = homotopHeur(SCNorm{sp,1}, SCMetr{sp,1}, avg_dist, heuristics, hScale);
    
    for i = 1:length(heuristics)                                            % iterate over types of heuristics
        for j = 1:hScale                                                    % iterate over homotopic scaling factors
            path = strcat('Results\sp',num2str(sparse(sp)), '\heur',num2str(heuristics(i)), '\h',num2str(j));
            mkdir(path);
                        
            hSC{i,j,sp} = normGraph(hSC{i,j,sp}, avg_roi_size, 'ROIsum', false, sparse(sp));
            hMetr{i,j,sp} = graphMetrics(hSC{i,j,sp}, path);
            
            % norm rows, i.e. input, to sum(CIJ) == 1:                      % remove this from connectivity data
            hSC{i,j,sp} = bsxfun(@rdivide,hSC{i,j,sp},sum(hSC{i,j,sp},2))'; % and make it part of the model?
            
            % loop over overall connection strength k, see Finger et al (2015) Fig.2 D
            for k = linspace(0.4,0.95,kScale);
                % RUN MODEL
                % acquire performances and store in #heuristics different h x k matrices
                
                % re-assess global and local prediction errors
                % get plots, organize results in reasonable folder structure
            end
        end
    end
    disp(strcat('job done: ',num2str(100*sp/length(sparse)),' %'))          % in case of impatience ...
end

save('Results\SC.mat','SCNorm', 'SCMetr');
save('Results\SCh.mat','hSC', 'hMetr');                                     % hSC is row-normalized    

% write a script for post-simulation analysis that plots metrics as fct of parameters