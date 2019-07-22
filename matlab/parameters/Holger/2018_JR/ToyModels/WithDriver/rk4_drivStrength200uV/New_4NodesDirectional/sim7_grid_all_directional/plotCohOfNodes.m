function plotCohOfNodes(FC_sum,node1,node2, d12, phase_offsets, path_results)
%PLOTCOHOFNODES Summary of this function goes here
%   Detailed explanation goes here

figure(1)
imagesc(d12,phase_offsets,squeeze(FC_sum(:,:,node1,node2))')
%set(gca,'clim',[0.3; 0.55]);
set(gca,'TickLength',[0 0]);
set(gca,'xTick',[0, 200]);
set(gca,'yTick',[0]);
xlabel('distance [mm]')
ylabel('stimulation phase offset [rad]')
set(gca,'YDir','normal')
title(['coh node ' num2str(node1) ' to node ' num2str(node2)])
colorbar;
saveas(gcf,fullfile(path_results,['coh_' num2str(node1) num2str(node2) '.png']))
saveas(gcf,fullfile(path_results,['coh_' num2str(node1) num2str(node2) '.pdf']))

end

