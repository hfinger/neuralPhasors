function patchCorrCov = fdr(patchCorrCov,p,alpha)

%alpha=0.01; %this is the FDR level

m=size(patchCorrCov,1)*size(patchCorrCov,2)*size(patchCorrCov,3);
c=sum(1./(1:m));

  pCurrent = p(:,:,:);
  patchCorrCovCurrent = patchCorrCov(:,:,:);
  [pCurrentSorted,pCurrentSortIds] = sort(pCurrent(:));
  tmp=c*m*pCurrentSorted./(1:m)';
  kplus1 = find(tmp > alpha, 1, 'first');
  patchCorrCovCurrent(pCurrentSortIds(kplus1:end)) = 0;
  patchCorrCov(:,:,:) = patchCorrCovCurrent;

