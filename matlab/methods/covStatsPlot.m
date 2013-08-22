function covStatsPlot( inStatsFile, rmax )
%FEATURECORRDENSMEANVSDISTEXPORT Summary of this function goes here
%   Detailed explanation goes here

featureStats=load(inStatsFile);

pos1 = [0.1 0.35 0.85 0.55];
pos2 = [0.1 0.05 0.74 0.25];
figure(1)
set(gcf,'PaperPositionMode', 'manual', ...
    'PaperUnits','centimeters', ...
    'Paperposition',[1 1 20 28.7])
  
  
  
  
subplot('Position',pos1)
imagesc(featureStats.meanStdOverAlpha(:,1:rmax)); 
colormap gray; 
colorbar; 
title('meanStdOverAlpha')
xlabel('radius [#pixel]'); 
ylabel('feature id');
subplot('Position',pos2)
plot(mean(featureStats.meanStdOverAlpha(:,1:rmax),1));
xlabel('radius [#pixel]'); 
ylabel('mean over Rfs');
print(gcf,'-dpdf','-r1200','meanStdOverAlpha.pdf'); 

subplot('Position',pos1)
imagesc(featureStats.meanStdOverRfs(:,1:rmax)); 
colormap gray; 
colorbar; 
title('meanStdOverRfs')
xlabel('radius [#pixel]'); 
ylabel('feature id');
subplot('Position',pos2)
plot(mean(featureStats.meanStdOverRfs(:,1:rmax),1));
xlabel('radius [#pixel]'); 
ylabel('mean over Rfs');
print(gcf,'-dpdf','-r1200','meanStdOverRfs.pdf'); 

if isfield(featureStats,'isotropy')
  subplot('Position',pos1) 
  imagesc(featureStats.isotropy(:,1:rmax)); 
  colormap gray; 
  colorbar; 
  title('isotropy')
  xlabel('radius [#pixel]'); 
  ylabel('feature id');
  subplot('Position',pos2) 
  plot(mean(featureStats.isotropy(:,1:rmax),1));
  xlabel('radius [#pixel]'); 
  ylabel('mean over Rfs');
  print(gcf,'-dpdf','-r1200','isotropy.pdf'); 
end

subplot('Position',pos1) 
imagesc(featureStats.meanAbsCorrDensVsRad(:,1:rmax)); 
colormap gray; 
colorbar; 
title('meanAbsCorr')
xlabel('radius [#pixel]'); 
ylabel('feature id');
subplot('Position',pos2) 
plot(mean(featureStats.meanAbsCorrDensVsRad(:,1:rmax),1));
xlabel('radius [#pixel]'); 
ylabel('mean over Rfs');
print(gcf,'-dpdf','-r1200','meanAbsCorrDensVsRad.pdf'); 

subplot('Position',pos1) 
imagesc(featureStats.meanCorrDensVsRad(:,1:rmax)); 
colormap gray; 
colorbar; 
title('meanCorr')
xlabel('radius [#pixel]'); 
ylabel('feature id');
subplot('Position',pos2) 
plot(mean(featureStats.meanCorrDensVsRad(:,1:rmax),1));
xlabel('radius [#pixel]'); 
ylabel('mean over Rfs');
print(gcf,'-dpdf','-r1200','meanCorrDensVsRad.pdf'); 

clf;
subplot(2,1,1) 
imagesc(featureStats.meanCorrDens);
colormap gray; 
colorbar; 
subplot(2,1,2) 
imagesc(featureStats.stdCorrDens);
colormap gray; 
colorbar; 
print(gcf,'-dpdf','-r1200','rfsCollapsedMeanCorrDens.pdf'); 

subplot(2,1,1) 
imagesc(featureStats.meanAbsCorrDens);
colormap gray; 
colorbar; 
subplot(2,1,2) 
imagesc(featureStats.stdAbsCorrDens);
colormap gray; 
colorbar; 
print(gcf,'-dpdf','-r1200','rfsCollapsedMeanAbsCorrDens.pdf'); 
% 
% clf;
% plot([featureStats.fractionCorrPos' featureStats.fractionCorrNeg']);
% xlabel('rf id')
% ylabel('fraction out/in connection probability')
% legend('excitatory','inhibitory')
% print(gcf,'-dpdf','-r1200','fractionConnProb.pdf'); 

end

