function plotGradientDescIter( savePath, pauseTime, iters, rfids, plotIds )
%PLOTGRADIENTDESCITER Summary of this function goes here
%   Detailed explanation goes here

if nargin<1 || isempty(savePath)
  savePath = pwd;
end
if nargin<2 || isempty(pauseTime)
  pauseTime=0.1;
end
if nargin<3 || isempty(iters)
  tmp=dir(fullfile(savePath,'gradientDescIter*.mat'));
  tmp=cellfun(@(x) str2double(x(17:end-4)), {tmp.name});
  iters = sort(tmp);
end
if nargin<5 || isempty(plotIds)
  plotIds = 1;
end

load(fullfile(savePath,'gradientDescStart.mat'),'featureOptions');
if featureOptions.usePCA
  preprocessing=load(fullfile(savePath,'preprocessing.mat'));
end

for iter=iters;
  load([savePath '/gradientDescIter' num2str(iter) '.mat']);
  if nargin<4 || isempty(rfids)
    rfids=1:size(W,1);
  end
  
  if featureOptions.biasUnit
    W = W(:,1:end-1);
  end
  
  if featureOptions.usePCA
    W = W*preprocessing.Vpca(1:featureOptions.PCAdim,:);
  end
  
  if any(plotIds==1)
    fh=figure(1);
    set(fh,'Name',['W Iter ' num2str(iter)]);
    [img] = plotColorFeatures( reshape(W',[featureOptions.patchsize featureOptions.patchsize 3 10 size(W,1)/10]), false, ['iter' num2str(iter) '.png'] );
    imshow(img,'InitialMagnification',300); title(['iter=' num2str(iter)]);
  end
  
  if any(plotIds==2)
    fh=figure(2);
    set(fh,'Name',['W Iter ' num2str(iter)]);
    [img] = plotColorFeatures( reshape(W',[featureOptions.patchsize featureOptions.patchsize 3 10 size(W,1)/10]), true, ['iter' num2str(iter) '.png'] );
    imshow(img,'InitialMagnification',300); title(['iter=' num2str(iter)]);
  end
  
  if any(plotIds==3)
    fh=figure(3);
    set(fh,'Name',['feature hist Iter ' num2str(iter)]);
    hist(reshape(feature(rfids,:),numel(feature(rfids,:)),1),-5:0.1:5);
    xlim([-5.1,5.1])
    
  end
  
  if any(plotIds==4)
    fh=figure(4);
    set(fh,'Name',['histW Iter ' num2str(iter)]);
    hist(reshape(W(rfids,:),numel(W(rfids,:)),1),-5:0.1:5);
    xlim([-6,6])
    
  end
  
  if any(plotIds==5)
    fh=figure(5);
    set(fh,'Name',['histA Iter ' num2str(iter)]);
    hist(reshape(A1(rfids,:),numel(A1(rfids,:)),1),0:0.1:1);
    xlim([-0.1,1.1])
  end
  
  if any(plotIds==6)
    colors=lines(5);
    fh=figure(6);
    set(fh,'Name',['Psi Iter ' num2str(iter)]);
    plot(PsiBefore(1:iter),'Color',colors(1,:))
    hold on;
    plot(PsiSparsetBefore(1:iter),'Color',colors(2,:))
    plot(PsiDecorrBefore(1:iter),'Color',colors(5,:))
    if exist('PsiMeanActBefore','var')
      plot(PsiMeanActBefore(1:iter),'Color',colors(3,:))
    end
    if exist('PsiSparserBefore','var')
      plot(PsiSparserBefore(1:iter),'Color',colors(4,:))
    end
    hold off;
    legend('Psi','PsiSparset','PsiDecorr','PsiMeanAct','PsiSparserBefore')
  end
  
    if any(plotIds==7)
      fh=figure(7);
      addC = featureOptions.sparserAddC;
      PsiSparser = (addC + mean(A1,1).^2) ./ (addC + mean(A1.^2,1));
      hist(PsiSparser);
    end

  pause(pauseTime);
end

end

