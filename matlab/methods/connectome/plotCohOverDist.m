function [ Coh_ordered ] = plotCohOverDist( dataStruct, stimPosID, maxDist, coherence_measure, dependent_vars, useShortestPath )
%PLOTCOHOVERDIST Shows coherence (colour-coded) of driven node and other nodes based on
%their distance to the driven node. Data will be ordered along y-axis
%according to the vals of dependent_var
%   Input Parameters:
%       dataStruct - structure with 1 field per file
%       stimPosID - ID of position of the driven node in ConnMat (scalar)
%       maxDist - maximum steps to take away from driven node (scalar)
%       coherence_measure - SE for Shannon Entropy or Coherence (string)
%       dependent_vars - 1 or 2 variables to sort data after (tuple)
%       useShortestPath - if true, make steps based on weighted shortest
%       paths. Else just plot respective row of ConnMat (logical)
%   Returns:
%       Coh_ordered - matrix that is plotted

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract relevant variables from structure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames = fieldnames(dataStruct);
n = length(fnames);
data_tmp = dataStruct.(fnames{1});
n_nodes = length(data_tmp.(coherence_measure));
stimNodes_tmp = data_tmp.stimPos;
n_drivers = length(stimNodes_tmp);
if ~iscell(dependent_vars)
    dvs_tmp = dependent_vars;
    dependent_vars = cell({});
    dependent_vars{1} = dvs_tmp;
end
n_dvs = length(dependent_vars);
if n_dvs > 2
    error('This function can plot coherence over distance for maximally two additional parameters')
end
m = min(maxDist, n_nodes-n_drivers);

% load SC mat and resort IDs
path_SCmat = '/net/store/nbp/projects/phasesim/databases/avg_SC.mat';
path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
load(path_SCmat)
load(path_ResortIDs)
clear D

% resort SC mat
C(1:33,:) = C(resortIdsMarlene,:);
C(34:end,:) = C(resortIdsMarlene + 33,:);
C(:,1:33) = C(:,resortIdsMarlene);
C(:,34:end) = C(:,resortIdsMarlene + 33);

% add homotopic connections
C = C + 0.1 * diag(ones(size(C,1)/2,1),size(C,1)/2) + 0.1 * diag(ones(size(C,1)/2,1),-size(C,1)/2);
C = bsxfun(@rdivide,C,sum(C,2));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract dependent variables and coherence measurement %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Coh_stimPos = zeros(n,m);
dvs = zeros(n,n_dvs);
    
for f=1:n
    
    % get relevant struct and driver positions
    data = dataStruct.(fnames{f});
    stimNodes = dataStruct.(fnames{f}).stimPos;
    stimPos = stimNodes(stimPosID);
    
    % get dvs out of struct
    for i=1:n_dvs
        dvs(f,i) = data.(dependent_vars{i})(1);
    end
    
    % get coherence measurement
    Coh = data.(coherence_measure);
    
    % extract coherence over distance
    if ~useShortestPath
        Coh_stimPos(f,:) = Coh(stimPos, n_drivers+1:end);
    else
        Coh_stimPos(f,:) = CohOverDist(Coh(n_drivers+1:end,n_drivers+1:end), C, stimPos, stimNodes, m);
    end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reorder coherence measurements according to dvs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n_dvs
    % look at unique values of each dv
    dv_unique = unique(dvs(:,i));
    n_vals = length(dv_unique);
    ind = 1:n;
    Coh_ordered = zeros(n_vals,m);
    % find positions where dvs correspond to unique vals
    for j=1:n_vals
        ind_tmp = ind(dvs(:,i) == dv_unique(j));
        %Coh_ordered(j,:) = Coh_stimPos(ind_tmp(round(length(ind_tmp)/2)),:); % find better choice than to take middle val
        Coh_ordered(j,:) = mean(Coh_stimPos(ind_tmp,:));
    end
    
    % plotting
    figure()
    y = linspace(dv_unique(1),dv_unique(end),length(dv_unique));
    x = 1:m;
    imagesc(x,y,Coh_ordered,[0,1])
    set(gca,'YDir','normal')
    colormap(jet)
    shading(gca,'interp')
    colorbar()
    title('Coherence of driven node over distance')
    xlabel('Distance from driven node')
    ylabel(dependent_vars{i})
    savefig(strcat('CohOverDist2D_',dependent_vars{i}))
end

if n_dvs > 1
    
    dv1 = unique(dvs(:,1));
    dv2 = unique(dvs(:,2));
    m1 = length(dv1);
    m2 = length(dv2);
    indA = 1:m1;
    indB = 1:m2;
    Sync = zeros(m1,m,m2);

    for j=1:n
        ind1 = indA(dv1 == dvs(j,1));
        ind2 = indB(dv2 == dvs(j,2));
        Sync(ind1,:,ind2) = Coh_stimPos(j,:);
    end
   
   figure()
   colormap('jet')
   [X,Y,Z] = meshgrid(1:m,dv1,dv2);
   h = slice(X,Y,Z,Sync, [], [], dv2);
   set(h, 'EdgeColor','none','FaceColor','interp')
   alpha(0.1)
   colorbar()
   title('Coherence over Distance and 2 DVs')
   xlabel('Distance')
   ylabel(dependent_vars{1})
   zlabel(dependent_vars{2})
   ylim([dv1(1),dv1(end)])
   zlim([dv2(1),dv2(end)])
   
   hold on
   pltDelay = 2;
   for k=1:length(h)    
       set(h(k),'FaceAlpha','interp','AlphaData',get(h(k),'ZData') + 1)
       drawnow()
       savefig(strcat('CohOverDist3D_',num2str(k)))
       pause(pltDelay)
   end
   alpha(0.3)
end

end

