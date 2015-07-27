function    RCsignif = metricsRC_wu(SC, ci, rcPerm, nPerm)
% function metricsRC_wu, gets called in metricsGlobal_wu
% SC = mean structural connectivity matrix
% ci = set of ca 19 individual connectivity matrices
% rcperm = set of rich-club vectors per permutated matrix: rcDim x noPermMatrices
% nPerm = number of group permutations to be performed

% tests significance of rich-club via permutation testing
% p-value 0.05, one-sided

% for references, see:  
%                       van den Heuvel et al. (2010)                        http://www.jneurosci.org/content/30/47/15915.full
%                       Lynall et al. (2010)                                http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2914251/
%                       Bassett et al. (2011)                               http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2878961/

%%% do individual networks ci still have to be preprocessed? >> YES! They're not symmetric etc
%%% die rcPerm sind permutationen auf SC -- sollten allerdings permutationen auf ci sein?

rcZ = @(x) (x-nanmean(rcPerm,2))./nanstd(rcPerm')';                        % compute z-score per individual network x 

cEmpty = cellfun(@isempty,ci);
% nIndiv = length(ci) - sum(cEmpty);                                       % equals sum(~cEmpty);

rcIndiv = zeros(66, sum(~cEmpty));
zIndiv  = zeros(66, sum(~cEmpty));
k = 1;

for j = 1:length(ci)                                                       % compute z-value of each individual network                                                       
    if (cEmpty(j)) continue; end
    rcIndiv(:,k) = rich_club_wu(ci{j},66);                                 % ab column 10 groesstenteils naN !!! << DO PREPROCESSING
    zIndiv(:,k)  = rcZ(rcIndiv(:,k));
    k = k + 1;
end

expBetw = nanmean(zIndiv,2);                                               % expBetw = between-group difference sample

permBetw = zeros(66,nPerm);
population = horzcat(rcIndiv, rcPerm);
for j = 1:nPerm                                                            % do permutation of grouped rich-club coefficients:                               
    [g1,idx] = datasample(population,sum(~cEmpty), 2,'Replace', false);    % draw new group 1
    g2 = population(:,setdiff(1:size(population,2),idx));                  % obtain a group 2
    
    rcPerm = g2;                                                           % overwrites g2 in anonymous fct rcZ, else: pass as arg
    for jj = 1:size(g1,2)                                                  % reassess between-group differences i.e. z-values     
       zIndiv(:,jj)  = rcZ(g1(:,jj));                                      % dim 66 x 19
    end
        
    permBetw(:,j) = nanmean(zIndiv,2);                                     % permBetw = between-group difference null distribution
end

% compare expBetw and permBetw
% calculate p-values

% Finally, for the S, lambda, and gamma parameters, the observed patient  
% versus control between-group effect (step 2) was assigned a p value by  
% computing the proportion of the total number of 10,000 entries resulting 
% from the permutation that was greater than (or smaller than if the effect 
% was negative) the observed patient versus control group effect. 
% A significance threshold of α = 0.05 (uncorrected) was used

compare = repmat(expBetw,1,nPerm) < permBetw;
frac    = sum(compare,2)./nPerm;

testii = 41;
RCsignif = 42;