load('net/store/nbp/projects/phasesim/results/Holger/JR_2Driver_drivScale.mat')
drivScales = 0:0.003:0.027;
phaseOffsets = linspace(0,2*pi,17);
dataRepeated = [data data(:,1)];

%%
figure(1);
clf;
plot(phaseOffsets,dataRepeated', 'Linewidth', 2)
legend(cellfun(@(x) num2str(x),num2cell(drivScales), 'UniformOutput', false), 'Location', 'NorthEastOutside')
ylabel('Coherence')
xlabel('tACS phase offset')
set(gca,'xtick',linspace(0,2*pi,5))
set(gca,'xticklabels',{0, '0.5\pi', '\pi', '1.5\pi', '2\pi'})
xlim([0 2*pi])
set(gca,'fontsize',14)
set(gca,'linewidth',2)

%%
figure(2);
surf(phaseOffsets, drivScales, dataRepeated)
xlabel('tACS phase offset')
set(gca,'xtick',linspace(0,2*pi,5))
set(gca,'xticklabels',{0, '0.5\pi', '\pi', '1.5\pi', '2\pi'})
xlim([0 2*pi])
zlabel('Coherence')
ylabel('tACS strength')
ylim([0 0.027])
set(gca,'fontsize',14)
set(gca,'linewidth',2)

%%
figure(3)
h = meshz(phaseOffsets, drivScales, dataRepeated);
xlabel('tACS phase offset')
set(gca,'xtick',linspace(0,2*pi,5))
set(gca,'xticklabels',{0, '0.5\pi', '\pi', '1.5\pi', '2\pi'})
xlim([0 2*pi])
zlabel('Coherence')
ylabel('tACS strength')
ylim([0 0.027])
set(gca,'fontsize',14)
set(gca,'linewidth',2)
set(h, 'faceColor', 'flat')
set(h, 'edgeColor', 'k')
axis vis3d
az = 55;
el = 40;
view(az, el);
% for t=0:0.02:21
%     az = 45 + 45*sin(t);
%     view(az, el);
%     drawnow;
% end
