results_path = '/net/store/nbp/projects/phasesim/results/Holger';

%% load old sim:

%load('net/store/nbp/projects/phasesim/results/Holger/JR_2Driver_drivScale.mat')
%drivScales = 0:0.003:0.027;
%numPhaseOffsets = 16;

%% OR load new sim:
numPhaseOffsets = 32;
numOtherParam = 26;
drivScales = 0:0.001:0.025;
numSim = numOtherParam*numPhaseOffsets;
% simResult = cell(numSim, 1);
% drivPosCoh = zeros(numSim, 1);
% for i=1:numSim
%   tmp = load(['/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_POvsScale/JR_2Driver_POvsScale_' num2str(i) '.mat']);
%   simResult{i} = tmp.simResult;
%   drivPosCoh(i) = tmp.simResult.drivPosCoh{1};
% end
% data = reshape(drivPosCoh, [numPhaseOffsets, numOtherParam])';
% save('/net/store/nbp/projects/phasesim/results/Holger/JR_2Driver_drivScale_new.mat','data','simResult')
load('/net/store/nbp/projects/phasesim/results/Holger/JR_2Driver_drivScale_new.mat')

%%
phaseOffsets = linspace(0,2*pi,numPhaseOffsets+1);
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
ylim([0 max(drivScales)])
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
ylim([0 max(drivScales)])
set(gca,'fontsize',14)
set(gca,'linewidth',2)
set(h, 'faceColor', 'flat')
set(h, 'edgeColor', 'k')
set(gca,'clim',[0 1]);
zData = get(get(gca,'Children'),'ZData');
zData(1,:) = 0;
zData(end,:) = 0;
zData(:,1) = 0;
zData(:,end) = 0;
set(get(gca,'Children'),'ZData',zData)
% axis vis3d
az = 55;
el = 40;
view(az, el);
% for t=0:0.02:21
%     az = 45 + 45*sin(t);
%     view(az, el);
%     drawnow;
% end

%% in-phase in center:

phaseOffsetsFliped = linspace(-pi,pi,numPhaseOffsets+1);
dataFlip = [data(:,16:end) data(:,1:16)];

%%
figh = figure(4);
h = meshz(phaseOffsetsFliped, drivScales, dataFlip);
% set(h,'LineWidth', 2)
xlabel('$\Delta \varphi$','Interpreter','LaTex')
set(gca,'xtick',linspace(-pi,pi,5))
set(gca,'xticklabels',{'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
set(gca,'ytick',[0 0.01, 0.02])
set(gca,'ztick',[0 0.5, 1])
xlim([-pi pi])
zlabel('Coherence')
ylabel('tACS')
ylim([0 max(drivScales)])
set(gca,'fontsize',20)
set(gca,'linewidth',3)
set(h, 'faceColor', 'flat')
set(h, 'edgeColor', 'k')
set(gca,'clim',[0 1]);
colormap jet;
zData = get(get(gca,'Children'),'ZData');
zData(1,:) = 0;
zData(end,:) = 0;
zData(:,1) = 0;
zData(:,end) = 0;
set(get(gca,'Children'),'ZData',zData)
shading interp
set(h','edgecolor','k')
az = 45;
el = 45;
view(az, el);

%% create video:
doCreateVideo = true;
% axis vis3d
if doCreateVideo
  daObj=VideoWriter(fullfile(results_path,'coh_surf')); %for default video format. 
  daObj.FrameRate=25;
%   axis vis3d;
%   camva('manual');
%   camva(8);
  open(daObj);
end
maxAngle = 45;
numInterp = 90;
win = [0.75+cos(linspace(0,pi,numInterp))/4 0.5*ones(1,100) 0.25+cos(linspace(0,pi,numInterp))/4 0*ones(1,100)];
win = [win 1-win];
for t=1:480%length(win)
  az = 90-90*win(t);
  el = 90*win(t);
  if az==0
    az=0.001;
    el=89.999;
  end
  view(az, el);
  drawnow;
  if doCreateVideo
    writeVideo(daObj,getframe(figh)); %use figure, since axis changes size based on view
  end
end

if doCreateVideo
  close(daObj);
end

%%
interesting = [0 0.006 0.011 0.012  0.013 0.018 0.02 0.025];
interestingIds = interesting*1000+1;
figure(5);
clf;
hold on;
cmap = jet;
hsv=rgb2hsv(cmap);
cm_data=interp1(linspace(0,1,size(cmap,1)),hsv,linspace(0,1,length(interesting)));
cm_data=hsv2rgb(cm_data);
for i=length(interesting):-1:1
  plot(phaseOffsetsFliped,dataFlip(interestingIds(i),:)', 'Linewidth', 2,'Color', cm_data(i,:));
end
legend(cellfun(@(x) num2str(x,'%.f m'),num2cell(1000*drivScales(interestingIds(end:-1:1))), 'UniformOutput', false), 'Location', 'EastOutside')
ylabel('Coherence')
xlabel('$\Delta \varphi$','Interpreter','LaTex')
set(gca,'xtick',linspace(0,2*pi,5))
set(gca,'xtick',linspace(-pi,pi,5))
set(gca,'xticklabels',{'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
xlim([-pi pi])
set(gca,'fontsize',14)
set(gca,'linewidth',2)

%% 
interesting = [0 0.006 0.011 0.012  0.013 0.018 0.02 0.025];
interestingLabels = {'0, flat','6, peak in-phase','11, peak anti-phase','12, peak 1/3 cycle','13, peak in-phase','18, peak in-phase','20, saturated in-phase','25, flat'};
interestingIds = interesting*1000+1;
figure(5);
clf;
hold on;
cmap = jet;
hsv=rgb2hsv(cmap);
cm_data=interp1(linspace(0,1,size(cmap,1)),hsv,linspace(0,1,length(interesting)));
cm_data=hsv2rgb(cm_data);
for i=length(interesting):-1:1
  plot(phaseOffsetsFliped,dataFlip(interestingIds(i),:)', 'Linewidth', 2,'Color', cm_data(i,:));
end
legend(interestingLabels(end:-1:1), 'Location', 'EastOutside')
ylabel('Coherence')
xlabel('$\Delta \varphi$','Interpreter','LaTex')
set(gca,'xtick',linspace(0,2*pi,5))
set(gca,'xtick',linspace(-pi,pi,5))
set(gca,'xticklabels',{'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
xlim([-pi pi])
ylim([0 1])
set(gca,'fontsize',14)
set(gca,'linewidth',2)