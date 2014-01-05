function plotMultCompare( stats, dim, groupNames )
%PLOTMULTCOMPARE Summary of this function goes here
%   Detailed explanation goes here
[c,m] = multcompare(stats,'dimension',dim);
% export_fig('populationMarginalMeansMethod.pdf', '-transparent', '-r72');
% barweb(m(:,1),1.96*m(:,2),0.6,[],[],'metrics','marginal mean',[],[],groupNames,[],'axis')

[c0p05] = multcompare(stats,'dimension',dim,'alpha',0.05);
[c0p01] = multcompare(stats,'dimension',dim,'alpha',0.01);
[c0p001] = multcompare(stats,'dimension',dim,'alpha',0.001);

figure; clf;
bar(m(:,1))
hold on;
h=errorbar(1:size(m,1),m(:,1),1.96*m(:,2),'c'); 
set(h,'linestyle','none')
set(h,'linewidth',2)
set(gca,'xticklabel',groupNames);
ylabel('marginal means')

pLevels = nan(1,size(c,1));
sig0p05 = abs(sign(c0p05(:,3))-sign(c0p05(:,5)))==0;
pLevels(sig0p05) = 0.05;
sig0p01 = abs(sign(c0p01(:,3))-sign(c0p01(:,5)))==0;
pLevels(sig0p01) = 0.01;
sig0p001 = abs(sign(c0p001(:,3))-sign(c0p001(:,5)))==0;
pLevels(sig0p001) = 0.001;
sigstar(num2cell(c(:,1:2),2),pLevels)

end

