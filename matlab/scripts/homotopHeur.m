function   hSC = homotopHeur(SC, SCMetr, heuristics, hScale, graphH0)

% implement heuristics to set homotopic connection weights

% effect of different sparseness values
% analyse graph for respective modification
% run model for each heuristics and compare performances
% what is the difference in effect on small-worldness for set/add ALL vs set/add MASKED?
% does a favorable heuristics candidate generalize on SC of finer parcellations? on graph incl subcortical ROI?
% compare set connections across heuristics via e.g. correlation, such as: corr(hDist, hInp) = 0.8510

h1 = logical(diag(ones(length(SC)/2,1),+length(SC)/2));                     % get masks for LH2RH and RH2LH homotopic connections
h2 = logical(diag(ones(length(SC)/2,1),-length(SC)/2));

% PARAMETERS:   (1) network sparseness                                      CHECK
%               (2) overall connection scaling k                            CHECK
%               (3) homotop connection scaling h
%               (4) rich-club connectedness, number of modules


%hSC = cell(length(heuristics),hScale,1);                                    % [# heuristics x resolution h]

switch heuristics
  case 1  % (a) set all weights (uniformly), see Messe et al (2014)
    add = false;                                                    % resets distribution
    %a = 5e-5;    b = 5e-4;                                         % use informed interval values or ...
    a = max(mean(mean(SC)) - hScale*mean(std(SC)), 0);              % ... set statistically, but
    b =     mean(mean(SC)) + hScale*mean(std(SC));                  % -> sensitive to % sparseness!
    
    
    hConn1 = 1;
    hConn2 = 1;
    
  case 2  % (b) add to all weights (uniformly), see Finger et al (2015)
    add = true;                                                     % preserves distribution
    
    a = 0.0;                                                        % h = linspace(0.0,0.22,12); --> linspace(0.02,.22, 11)
    b = 0.22;                                                       % 'fraction of additional homotopic input per ROI'
    
    hInp1 = SCMetr.perROI.inpStrength(1:length(SC)/2);
    hInp2 = SCMetr.perROI.inpStrength(length(SC)/2 + 1:length(SC));
    hInp  = mean(horzcat(hInp1, hInp2), 2);                         % take mean to conserve symmetry
    
    hConn1 = hInp;                                                  % see ConnectomeSim.m, line 240
    hConn2 = hInp;
    
  case 3  % (c) set all weights (euclid. dist. between homologous regions)
    add = false;
    a = 0.0;   b = 0.22;
    
    hDist1 = 1./SCMetr.perConn.euclDist(h1);                        % map euclidean distance to connection weight,
    hDist2 = 1./SCMetr.perConn.euclDist(h2);                        % avg_dist is not symmetric, why? Pearson corr' .999
    hDist  = mean(horzcat(hDist1, hDist2), 2);
    
    hConn1 = hDist;                                                 % interestingly, corr(hDist, hInp) = 0.8510
    hConn2 = hDist;                                                 % (Pearson's ignores difference in magnitude, but: is comparable!)
    
  case 4  % (d) set all weights (participation index)
    add = false;
    a = 0.05;   b = 0.6;
    
  case 5  % (e) set weights between RC members (mask) -- requires > 0 permutation graphs
    if graphH0 > 0
      add = false;
      a = 0.05;   b = 0.6;
      
      rc = min(SCMetr.perGlob.Rsize(~mod(SCMetr.perGlob.Rsize,2)));   % choose a particular RC: smallest-size even numbered
      hRC = SCMetr.perGlob.rcMemberS(:,find(SCMetr.perGlob.Rsize==rc,1,'first'));
      hRC = reshape(hRC, length(SC)/2,2);
      hRC = or(hRC(:,1),hRC(:,2));                                    % apply OR or AND ?
      
      hConn1 = hRC;
      hConn2 = hRC;
    end
    
  case 6  % (f) set weights between connector hubs (mask)
    add = false;
    a = 0.05;   b = 0.6;
    
  case 7  % (g) compute participation index (or other metric) for each ROI in
    % ipsilateral hemisphere, then compare values with homologous ROI across hemisphere --
    % P1: structural embedding of ROI determines its function, Passingham et al (2002) 'Anat. Basis Funct. Localization'
    % P2: brain asymmetry / lateralization := differences in function
    % C : ...
    % differences in value might hint at differences in function,
    % --> brain asymmetry / lateralized function
    add = false;
    a = 0.00;   b = 0.022;                                              % metrics here are normalized/within [0,1]
    % b = polyfit(MIdSimi(h1),0.22*hInp,0)
    
    hLeft = SC(1:length(SC)/2,1:length(SC)/2);
    hRight= SC(length(SC)/2+1:length(SC),length(SC)/2+1:length(SC));
    
    % http://brenocon.com/blog/2012/03/cosine-similarity-pearson-correlation-and-ols-coefficients/
    cosSimi = (dot(hLeft,hRight) ./ (norm(hLeft,2) *norm(hRight,2)))';  % cosine similarity [0,1] for SC>0
    jacSimi = (sum(min(hLeft,hRight),1) ./ sum(max(hLeft,hRight),1))';  % generalized Jaccard similarity
    corSimi = corr(hLeft, hRight);                                      % Pearson's linear correlation coeff
    corSimi = corSimi(eye(length(SC)/2)==1);                            % ALL ROI highly +correlated, useless?
    
    MIdSimi = matching_ind_und(SC>0);                                   % matching index (on binary graph), req. sparseness
    % MId is symmetric, hence MIdSimi(h1) = MIdSimi(h2)
    
    hConn1 = MIdSimi(h1);
    hConn2 = MIdSimi(h2);
    
    % metacorr: 1) corr(1./avg_dist(h1), MIdSimi(h1)) = 0.55
    %           2) corr(cosSimi,jacSimi) = 0.37
    %           3) corr(1./avg_dist(h1), corSimi) = 0.29
    %           4) corr(1./avg_dist(h1), cosSimi) = 0.85
    
  otherwise % do not change homotopic connections (default SC)
    add = true;
    a = 0;   b = 0;
    
    hConn1 = 1;
    hConn2 = 1;
end

% if no permutation graphs provided, set heuristics parameters to generic values
if ~(graphH0>0)
  add = false;
  a = 0.0;      b = 0.22;
  hConn1 = 1;   hConn2 = 1;
end

% perform ADD or SET of homotopic connections, then scale ALL or MASKED
maxRes  = 12;
scaling = linspace(a,b,maxRes);

SC(h1) = SC(h1)*add + scaling(min(12,hScale)) * hConn1;
SC(h2) = SC(h2)*add + scaling(min(12,hScale)) * hConn2;

% inserting homotopic connections affects sparseness! how to control for this?
effSparse = sum(sum(SC == 0)) / numel(SC);                          % determine effective sparseness of resulting graph

hSC = SC;

% issues:   - introducing homotopic connections affects sparseness
%             --> this is naturally going to affect graph metrics, CONTROL FOR THIS! (random permutation?)
%           - how to choose reasonable scaling intervals h for fair comparison? -> b=polyfit(hConn1,0.22*hInp,0)