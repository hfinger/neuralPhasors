function covStats( inCorrFile, outStatsFile, calcIsotropy, calcIsotropyCov, filterFcn, FDRalpha )
%FEATURECORRDENSMEANVSDIST Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
  calcIsotropy = false;
end
if nargin<4
  calcIsotropyCov = false;
end

corrs = load(inCorrFile);

numFeatures = size(corrs.patchCorrCov,3);
rmax = (size(corrs.patchCorrCov,1)-1)/2;

dx=repmat((-rmax:rmax)',[1 2*rmax+1]);
dy=repmat(-rmax:rmax,[2*rmax+1, 1]);
r = floor(0.5+sqrt(dx.^2+dy.^2));
r=r(:);
ids = cell(rmax,1);
for rr=1:rmax
  ids{rr} = find(r==rr);
end

meanCorrDens = sum(sum(corrs.patchCorrCov,3),4) / numFeatures^2;
meanCorrDensSqr = sum(sum(corrs.patchCorrCov.^2,3),4) / numFeatures^2;
meanAbsCorrDens = sum(sum(abs(corrs.patchCorrCov).^2,3),4) / numFeatures^2;
fractionCorr = squeeze(sum(sum(corrs.patchCorrCov,1),2));
fractionCorr = bsxfun(@rdivide,fractionCorr,sum(fractionCorr,2));

isotropy=zeros(numFeatures,rmax);
isotropyCov=zeros(numFeatures,rmax);
meanStdOverAlpha=zeros(numFeatures,rmax);
meanCoefOfVariationOverAlpha=zeros(numFeatures,rmax);
meanStdNormalizedOverAlpha=zeros(numFeatures,rmax);
meanStdOverRfs=zeros(numFeatures,rmax);
meanCorrDensVsRad=zeros(numFeatures,rmax);
meanAbsCorrDensVsRad=zeros(numFeatures,rmax);

oneBatch.meanStdOverAlpha=cell(1,rmax);
oneBatch.meanCoefOfVariationOverAlpha=cell(1,rmax);
oneBatch.meanStdNormalizedOverAlpha=cell(1,rmax);
oneBatch.meanStdOverRfs=cell(1,rmax);
oneBatch.meanCorrDensVsRad=cell(1,rmax);
oneBatch.meanAbsCorrDensVsRad=cell(1,rmax);

for f0id=1:numFeatures
  disp(num2str(f0id))
  
  if nargin>5 && ~isempty(FDRalpha)
    corrs.patchCorrCov(:,:,:,f0id) = fdr(corrs.patchCorrCov(:,:,:,f0id),corrs.p(:,:,:,f0id),FDRalpha);
  end
  
  %% calculate isotropy index of feature
  corrDens = reshape(corrs.patchCorrCov(:,:,:,f0id),[size(corrs.patchCorrCov,1)*size(corrs.patchCorrCov,2) size(corrs.patchCorrCov,3)]);
  
  if nargin>4 && ~isempty(filterFcn)
    corrDens = feval(filterFcn,corrDens);
  end
  
  for rr=1:rmax
    if calcIsotropy
      RHO = corr(corrDens(ids{rr},:)');
      RHO = RHO(:);
      RHO(1:(length(ids{rr})+1):length(ids{rr})^2)=[];
      isotropy(f0id,rr)=mean(RHO);
    end
    
    if calcIsotropyCov
      RHO = cov(corrDens(ids{rr},:)');
      RHO = RHO(:);
      RHO(1:(length(ids{rr})+1):length(ids{rr})^2)=[];
      isotropyCov(f0id,rr)=mean(RHO);
    end
    
    tmp = corrDens(ids{rr},:);
    
    
    meanStdOverAlpha(f0id,rr)=mean(std(tmp,[],1));
    meanCoefOfVariationOverAlpha(f0id,rr)=mean(std(tmp,[],1) ./ mean(tmp,1));
    meanStdNormalizedOverAlpha(f0id,rr)=mean(std(tmp,[],1) ./ mean(abs(tmp),1));
    meanStdOverRfs(f0id,rr)=mean(std(tmp,[],2));
    meanCorrDensVsRad(f0id,rr)=mean(tmp(:));
    meanAbsCorrDensVsRad(f0id,rr)=mean(abs(tmp(:)));
    
    oneBatch.meanStdOverAlpha{rr}{f0id}=(std(tmp,[],1));
    oneBatch.meanCoefOfVariationOverAlpha{rr}{f0id}=(std(tmp,[],1) ./ mean(tmp,1));
    oneBatch.meanCoefOfVariationOverAlpha{rr}{f0id}(isnan(oneBatch.meanCoefOfVariationOverAlpha{rr}{f0id})) = 0;
    oneBatch.meanStdNormalizedOverAlpha{rr}{f0id}=(std(tmp,[],1) ./ mean(abs(tmp),1));
    oneBatch.meanStdNormalizedOverAlpha{rr}{f0id}(isnan(oneBatch.meanStdNormalizedOverAlpha{rr}{f0id})) = 0;
    oneBatch.meanStdOverRfs{rr}{f0id}=(std(tmp,[],2))';
    oneBatch.meanCorrDensVsRad{rr}{f0id}=(tmp(:))';
    oneBatch.meanAbsCorrDensVsRad{rr}{f0id}=(abs(tmp(:)))';
  end
  
end

oneBatch.meanStdOverAlpha = cell2mat(cellfun(@(x) nanmean(cell2mat(x)), oneBatch.meanStdOverAlpha,'UniformOutput',false));
oneBatch.meanCoefOfVariationOverAlpha = cell2mat(cellfun(@(x) nanmean(cell2mat(x)), oneBatch.meanCoefOfVariationOverAlpha,'UniformOutput',false));
oneBatch.medianStdNormalizedOverAlpha = cell2mat(cellfun(@(x) median(cell2mat(x)), oneBatch.meanStdNormalizedOverAlpha,'UniformOutput',false));
oneBatch.meanStdNormalizedOverAlpha = cell2mat(cellfun(@(x) nanmean(cell2mat(x)), oneBatch.meanStdNormalizedOverAlpha,'UniformOutput',false));
oneBatch.meanStdOverRfs = cell2mat(cellfun(@(x) nanmean(cell2mat(x)), oneBatch.meanStdOverRfs,'UniformOutput',false));
oneBatch.meanCorrDensVsRad = cell2mat(cellfun(@(x) nanmean(cell2mat(x)), oneBatch.meanCorrDensVsRad,'UniformOutput',false));
oneBatch.meanAbsCorrDensVsRad = cell2mat(cellfun(@(x) nanmean(cell2mat(x)), oneBatch.meanAbsCorrDensVsRad,'UniformOutput',false));

twoMeans.meanStdOverAlpha = nanmean(meanStdOverAlpha,1);
twoMeans.meanCoefOfVariationOverAlpha = nanmean(meanCoefOfVariationOverAlpha,1);
twoMeans.medianStdNormalizedOverAlpha = median(meanStdNormalizedOverAlpha,1);
twoMeans.meanStdNormalizedOverAlpha = nanmean(meanStdNormalizedOverAlpha,1);
twoMeans.meanStdOverRfs = nanmean(meanStdOverRfs,1);
twoMeans.meanCorrDensVsRad = nanmean(meanCorrDensVsRad,1);
twoMeans.meanAbsCorrDensVsRad = nanmean(meanAbsCorrDensVsRad,1);

varCorrDens = meanCorrDensSqr - meanCorrDens.^2;
varAbsCorrDens = meanCorrDensSqr - meanAbsCorrDens.^2;
stdCorrDens = sqrt(varCorrDens);
stdAbsCorrDens = sqrt(varAbsCorrDens);
fractionCorr = sum(fractionCorr,1);

save(outStatsFile,'oneBatch','twoMeans','meanCorrDens','meanCorrDensSqr','varCorrDens','stdCorrDens','meanAbsCorrDens','varAbsCorrDens','stdAbsCorrDens','meanStdOverAlpha','meanCoefOfVariationOverAlpha','meanStdNormalizedOverAlpha','meanStdOverRfs','meanCorrDensVsRad','meanAbsCorrDensVsRad','fractionCorr');
if calcIsotropy
  save(outStatsFile,'-append','isotropy');
end
if calcIsotropyCov
  save(outStatsFile,'-append','isotropyCov');
end


end

