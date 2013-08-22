Y=simResult.Y;
numCluster=network.numCluster;
N=network.N;
f=sim.f;
v=sim.v;

%%
delayInRad = 2*pi*f*D/(v*1e3);
pariedSync = bsxfun( @minus, reshape(Y,[1 size(Y,1) size(Y,2)]), reshape(Y,[size(Y,1) 1 size(Y,2)]) );
pariedSync = cos(bsxfun( @minus, pariedSync, delayInRad));

figure(7);
for t=1:10:length(localPhase)
  clf;
  
  subplot(1,2,1);
  plotPolarPhases( localPhase(:,t), numCluster, N );
  title(['time: ' num2str(t)])
%   drawnow;
  
  subplot(1,2,2);
  currSync = pariedSync(:,:,t);
  currSync(C==0) = NaN;
  %   imagesc(currSync)
  colormap jet;
  pcolor([currSync nan(size(currSync,2),1); nan(1,size(currSync,1)+1)]);
  colormap jet;
  colorbar;
  hold all; 
  shading flat;
  axis square
  set(gca, 'ydir', 'reverse');
  xlabel('source neuron')
  ylabel('target neuron')
  set(gca,'clim',[-1 1]);
  title(['time: ' num2str(t)])
  drawnow;
  
  
  
  
end