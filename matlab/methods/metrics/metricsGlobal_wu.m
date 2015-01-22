function   perGlob = metricsGlobal_wu(SC, nSamples)

paths = dataPaths();
load(fullfile(paths.databases,'SC_Bastian','Stuct_connectivity_for_Holger.mat'), 'struct_labels');

w2d = @(M) (1./M);                                                          % define mapping from weight matrix to distance-length matrix
% w2d = @(M) (-log(M));                                                     % alternative mapping, applied in Goni et al (2014)

perGlob.clCoeff = mean(clustering_coef_wu(SC));
[perGlob.lambda, perGlob.efficiency] = charpath(distance_wei(w2d(SC)));     % avg shortest path length,
% avg inv shortest path length, same as efficiency_wei(CIJ, 0)

[perGlob.R] = rich_club_wu(SC, 66);                                         % Opsahl et al. Phys Rev Lett, 2008, 101(16)
%%%%%%[CIJscore,sn] = score_wu(SC,s);                                       % s-core decomposition, see van den Heuvel & Sporns (2011), Fig. 5

perGlob.sparse = sum(sum(SC == 0)) / numel(SC);                             % determine level of sparseness, [0, 1]

%% get normalized measurements

perm = zeros(nSamples, 2);
rcPerm = zeros(66, nSamples);
for j = 1:nSamples                                                          % generate a total of nSamples H0 networks:
  %W0 = randmio_und_connected(SC, 5);                                       % see Maslov & Sneppen (2002) -- preserves node degrees
  [W0,~] = null_model_und_sign(SC);                                         % see Rubinov & Sporns (2011) -- preserves node degrees + weights
  perm(j,1) = mean(clustering_coef_wu(W0));                                 % CC IS ALWAYS THE SAME for Rubinov & Sporns -- WHY?
  perm(j,2) = charpath(distance_wei(w2d(W0)));
  rcPerm(:,j) = rich_club_wu(W0,66);
end

% compute normalized graph metrics;
% compare to van den Heuvel & Sporns (2011), Table 1

randClus = mean(perm(:,1),1);
perGlob.normC = perGlob.clCoeff/randClus;                                   % get normalized clustering coefficients

randChar = mean(perm(:,2),1);
perGlob.normP = perGlob.lambda/randChar;                                    % get normalized global efficiency

% Sporns et al (2004) summarized 'structural cortical networks have revealed small-world attributes, with path lengths
% that are close to those of equivalent random networks but with significantly higher values for the clustering coefficient'

% Cexp = 9e-5, Crand = 6e-5
% Eexp = 5e-4, Erand = 6e+3

perGlob.small = perGlob.normC/perGlob.normP;                                % get small-world ratio, see Humphries & Gurney (2008)

randRich = nanmean(rcPerm,2);
perGlob.Rnorm = perGlob.R ./ randRich';                                     % get normalized rich-club coefficient

% perform 1-tailed test (parameter-free) for significance of richt-club topology via permutation testing

% if nargin == 3
%     perGlob.RC = metricsRC_wu(SC, ci, rcPerm, 1000);                      % literature samples 5k-10k permutations to get H0 distribution
% end

alpha   = 0.05;                                                             % van den Heuvel et al (2010) used alpha = 0.05 (uncorrected)
% van den Heuvel & Sporns (2011) used alpha = 0.001 (subcort ROI incl)
compare = repmat(perGlob.R',1,nSamples) > rcPerm;                           % use perGlob.Rnorm instead !
pVal    = 1-sum(compare,2)./nSamples;
perGlob.Rsignif  = (pVal < alpha)';

degrees = sum((SC~=0));
rcMember = bsxfun(@ge,degrees',(1:length(degrees)));
rcSize   = sum(rcMember,1);

perGlob.Rsize = rcSize(perGlob.Rsignif);                                    % get sizes of significant RC
perGlob.rcMemberS = rcMember(:,perGlob.Rsignif');                           % get members of significant RC
%struct_labels=repmat(struct_labels,1,sum(perGlob.Rsignif));

for j = 1:sum(perGlob.Rsignif)
  perGlob.Rlabels{:,j} = struct_labels(perGlob.rcMemberS(:,j));             % get labels of significant RC
end

% check if RC membership is symmetric to allow for straightforward homotopic connections

dbpoint = 42;
