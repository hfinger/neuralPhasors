function   hSC = homotopHeur(SC, SCMetr, heuristics, hScale, graphH0)

% implement heuristics to set homotopic connection weights

% effect of different sparseness values
% analyse graph for respective modification
% run model for each heuristics and compare performances
% what is the difference in effect on small-worldness for set/add ALL vs set/add MASKED?
% does a favorable heuristics candidate generalize on SC of finer parcellations? on graph incl subcortical ROI?
% compare set connections across heuristics via e.g. correlation, such as: corr(hDist, hInp) = 0.8510

h1 = logical(diag(ones(length(SC)/2,1),+length(SC)/2));                     % get masks for LH2RH and
h2 = logical(diag(ones(length(SC)/2,1),-length(SC)/2));                     % RH2LH homotopic connections

% PARAMETERS:   (1) network sparseness                                      CHECK
%               (2) overall connection scaling k                            CHECK
%               (3) homotop connection scaling h                            individual per heuristics
%                   ---> validity of individual scaling can be checked by HvsK performance plot, local maximum? blob?
%               (4) rich-club connectedness, number of modules

switch heuristics
  case 1  % (a) set all weights (uniformly), see Messe et al (2014)
    add = false;                                                    % resets distribution
    %a = 5e-5;    b = 5e-4;                                         % use informed interval values or ...
    a = max(mean(mean(SC)) - 2*mean(std(SC)), 0);                   % ... set statistically, but
    b =     mean(mean(SC)) + 2*mean(std(SC));                       % -> sensitive to % sparseness and size of graph!
   
    hConn1 = 1;                                                     % hConn = 1 if add = true
    hConn2 = 1;                                                     %       = ... else
    
  case 2  % (b) add to all weights (uniformly), see Finger et al (2015)
    add = true;                                                     % preserves distribution
    
    a = 0.0; % 0.02;                                                % h = linspace(0.0,0.22,12); --> linspace(0.02,.22, 11)
    b = 0.1; % 0.22;                                                % 'fraction of additional homotopic input per ROI'
    
    % SC = bsxfun(@rdivide,SC,SCMetr.perROI.inpStrength);           % IF b=.22, THEN normalize AND set hConn = 1
    
    hInp1 = SCMetr.perROI.inpStrength(1:length(SC)/2);
    hInp2 = SCMetr.perROI.inpStrength(length(SC)/2+1:length(SC));
    hInp  = mean(horzcat(hInp1, hInp2), 2);                         % take mean to conserve symmetry
    
    hConn1 = hInp;                                                  % see ConnectomeSim.m, line 240
    hConn2 = hInp;
    
  case 3  % (c) set all weights (euclid. dist. between homologous regions)
    add = true;
    a = 0.0;   b = 1e-4;                                            % b=1e-3 for add=false, b=1e-4 for add=true
                                                                    % perf .67 for add=false, for add=true
    
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
    add = true;
    a = 0.00;   b = 1e-4;                                               % metrics here are normalized/within [0,1]
    % b = polyfit(MIdSimi(h1),0.22*hInp,0)                              % b=1e-5 for add=true
    
    hLeft = SC(1:length(SC)/2,1:length(SC)/2);
    hRight= SC(length(SC)/2+1:length(SC),length(SC)/2+1:length(SC));
    
    hConn1 = -log(SCMetr.perROI.cosSimiC); %cosSimiR(h1);
    hConn2 = -log(SCMetr.perROI.cosSimiC); %cosSimiR(h2);
    
    % metacorr: 1) corr(1./avg_dist(h1), MIdSimi(h1)) = 0.55
    %           2) corr(cosSimi,jacSimi) = 0.37
    %           3) corr(1./avg_dist(h1), corSimi) = 0.29
    %           4) corr(1./avg_dist(h1), cosSimi) = 0.85
    
    % cosSimi: b = 1e-4;
    % MIdSimi: b = 5e-6;
    
    % perf_MIdSimi compares to perf_h2
    
  otherwise % do not change homotopic connections (default SC)
    add = true;
    a = 0;   b = 0;
    
    hConn1 = 1;
    hConn2 = 1;
end

% if no permutation graphs provided, set heuristics parameters to generic values
if ~(graphH0>0) && heuristics == (5|6)
  add = false;
  a = 0.0;      b = 0.22;
  hConn1 = 1;   hConn2 = 1;
end

% perform ADD or SET of homotopic connections, then scale ALL or MASKED
maxRes  = 12;
scaling = linspace(a,b,maxRes);

SC(h1) = SC(h1)*add + scaling(min(maxRes,hScale)) * hConn1;
SC(h2) = SC(h2)*add + scaling(min(maxRes,hScale)) * hConn2;

% inserting homotopic connections affects sparseness! how to control for this?
effSparse = sum(sum(SC == 0)) / numel(SC);                          % determine effective sparseness of resulting graph

hSC = SC;

% issues:   - introducing homotopic connections affects sparseness
%             --> this is naturally going to affect graph metrics, CONTROL FOR THIS! (random permutation?)
%           - how to choose reasonable scaling intervals h for fair comparison? -> b=polyfit(hConn1,0.22*hInp,0)