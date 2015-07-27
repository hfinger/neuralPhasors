function   perGlob = metricsGlobal_wd(CIJ, nSamples)

perGlob.clCoeff = mean(clustering_coef_wd(CIJ));                           % see Fagiolo (2007) Phys Rev E 76:026107
[perGlob.lambda, perGlob.efficiency] = charpath(distance_wei(1./CIJ));     % avg shortest path length,
% avg inv shortest path length, same as efficiency_wei(CIJ, 0)

[perGlob.R] = rich_club_wd(CIJ,66);                                        % Opsahl et al. Phys Rev Lett, 2008, 101(16)

%% get normalized measurements

if nSamples >0
  perm = zeros(nSamples, 2);
  rcperm = zeros(66, nSamples);
  for j = 1:nSamples
    [W0,~] = null_model_dir_sign(CIJ);                                     % generate null-hypothesis network, see Rubinov & Sporns (2011)
    perm(j,1) = mean(clustering_coef_wd(W0));                              % CC IS ALWAYS THE SAME -- WHY?
    perm(j,2) = charpath(distance_wei(1./W0));
    rcperm(:,j) = rich_club_wd(W0,66);
  end
  
  % compute normalized graph metrics;
  % compare to van den Heuvel & Sporns (2011), Table 1
  
  randClus = mean(perm(:,1),1);
  perGlob.normC = perGlob.clCoeff/randClus;                                  % get normalized clustering coefficients
  
  randChar = mean(perm(:,2),1);
  perGlob.normP = perGlob.lambda/randChar;                                   % get normalized global efficiency
  
  perGlob.small = perGlob.normC/perGlob.normP;                               % get small-world ratio, see Humphries & Gurney (2008)
  
  randRich = mean(rcperm,2);
  perGlob.normR = perGlob.R ./ randRich';                                    % get normalized rich-club coefficient
end