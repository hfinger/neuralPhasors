function plotCovAndConn( inCorrConnFile, inRfWeightsFile, featureIds, outFilePrefix, maxLength, padLength )
%PLOTCOVANDCONN Summary of this function goes here
%   Detailed explanation goes here

corrs = load(inCorrConnFile);
corrs2 = zeros(2*maxLength+1, 2*maxLength+1, length(featureIds), length(featureIds));

if isfield(corrs,'W')
  ids = ismember( corrs.W.f0, featureIds) & ismember( corrs.W.f1, featureIds);
  f0 = corrs.W.f0(ids);
  f1 = corrs.W.f1(ids);
  dx = corrs.W.dx(ids);
  dy = corrs.W.dy(ids);
  w = corrs.W.w(ids);
  for connid = 1:length(f0)
    xind = maxLength+1-dx(connid);
    yind = maxLength+1-dy(connid);
    if xind>0 && xind<=size(corrs2,1) && yind>0 && yind<=size(corrs2,2)
      corrs2(xind,yind,f0(connid)==featureIds,f1(connid)==featureIds) = w(connid);
    end
  end
elseif isfield(corrs,'patchCorrCov')
  sizeDiff = (size(corrs.patchCorrCov,1) - size(corrs2,1))/2;
  corrs2 = corrs.patchCorrCov(sizeDiff+1:end-sizeDiff,sizeDiff+1:end-sizeDiff,featureIds,featureIds);
end
plotFeatures( corrs2 , [outFilePrefix 'Corrs.tikz'], 'gray', padLength, [],[],[],[],[],true);
axis off;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
print(gcf,'-dpng',[outFilePrefix '.png']);
print(gcf,'-dpsc2',[outFilePrefix '.eps']);
% matlab2tikz( [outFilePrefix '.tikz'] , 'width', '\figurewidth', 'relativePngPath', '', 'interpretTickLabelsAsTex', true );

if ~isempty(inRfWeightsFile)
  weights = load(inRfWeightsFile);
  % patches = zeros(2*maxLength+1, 2*maxLength+1, 3, length(featureIds));
  % sizeDiff1 = ceil((size(patches,1) - size(weights.W,1))/2);
  % sizeDiff2 = floor((size(patches,1) - size(weights.W,1))/2);
  % patches(sizeDiff1+1:end-sizeDiff2,sizeDiff1+1:end-sizeDiff2,:,:,:) = squeeze(weights.W(:,:,:,1,1,featureIds));
  % imwrite(patches,[outFilePrefix 'Horz.png'],'png')

  patches = permute(weights.W(:,:,:,1,1,featureIds),[1 2 3 6 4 5]);
  % patches = cat(5,patches,zeros(size(patches)));
  plotColorFeatures( patches , true, [outFilePrefix 'Horz.png'], true, (size(corrs2,1)+padLength)-size(patches,1) );

  patches = permute(weights.W(:,:,:,1,1,featureIds),[1 2 3 4 6 5]);
  % patches = cat(4,patches,zeros(size(patches)));
  plotColorFeatures( patches , true, [outFilePrefix 'Vert.png'], true, (size(corrs2,1)+padLength)-size(patches,1) );
end

% if ~isempty(inConnWeightsFile)
%   conn = load(inConnWeightsFile);
%   connPattern = zeros(size(corrs2));
%   syncIds=find(conn.W.w==1);
%   desyncIds=find(conn.W.w==-1);
%   dxId = conn.W.dx+(size(corrs2,1)-1)/2;
%   dyId = conn.W.dy+(size(corrs2,2)-1)/2;
%   connPattern(sub2ind(size(corrs2), dxId(syncIds), dyId(syncIds), conn.W.f0(syncIds), conn.W.f1(syncIds))) = 1;
%   connPattern(sub2ind(size(corrs2), dxId(desyncIds), dyId(desyncIds), conn.W.f0(desyncIds), conn.W.f1(desyncIds))) = -1;
%   
%   plotFeatures( connPattern , [outFilePrefix 'Conns.tikz'], 'gray', padLength);
%   axis off;
%   set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
%   print(gcf,'-dpng',[outFilePrefix 'Conns.png']);
%   print(gcf,'-dpsc2',[outFilePrefix 'Conns.eps']);
% end

end

