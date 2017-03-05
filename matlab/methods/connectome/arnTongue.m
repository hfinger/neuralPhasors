function [ sync ] = arnTongue( dataStruct, targets, coherence_measure, cohInd, dependent_vars, clusts )
% ARNTONGUE function that loops through all entries in dataStruct and plots arnold
% tongue between 2 target regions, i.e. their respective coherence measure
% in dependence of dependent_vars
%   Input Parameters:
%       dataStruct - structure with 1 entry for each file
%       targets - driving and driven region in connectivity matrix (tuple)
%       coherence_measure - SE for Shannon Entropy or Coherence (string)
%       dependent_vars - 2 dvs for which to plot AT (cell of strings)
%       clusts - if Connectivity matrix is clustered, enter cluster
%       indices, else enter 0
%   Returns:
%       sync - array which is plotted as AT

%%
% find fieldnames and initialize stuff
fnames = fieldnames(dataStruct);
n = length(fnames);
dvs = zeros(n,2);

%%
% get simulation parameters that are the same for all files
data_tmp = dataStruct.(fnames{1});
C = data_tmp.(coherence_measure){1,1};
n_nodes = length(C);

%%
% make clusters
if clusts == 0
    clusts = cell(1,n_nodes);
    for i = 1:n_nodes
        clusts{i} = i;
    end
end

%%
% loop through data files and evaluate synchronization
S = zeros(n, length(clusts), length(clusts));
for f = 1:n
    
    data = dataStruct.(fnames{f});
    Sync = squeeze(data.(coherence_measure){1,cohInd});

    
    % get dependent variables from struct
    dv1 = data.(dependent_vars{1});
    if length(dv1) > 1
        dvs(f,1) = data.(dependent_vars{1})(targets(1));
    else
        dvs(f,1) = data.(dependent_vars{1});
    end
    
    dv2 = data.(dependent_vars{2});
    if length(dv2) > 1
        dvs(f,2) = data.(dependent_vars{2})(targets(1));
    else
        dvs(f,2) = data.(dependent_vars{2});
    end
    
    % evaluate synchronization between clusters
    for r = 1:length(clusts)
        for c = 1:length(clusts)
            S(f,r,c) = mean(reshape(Sync(clusts{r},clusts{c}),[1,length(clusts{r}) * length(clusts{c})]));
        end
    end
    
end

%%
% sort synchronization data according to stim freqs and strengths and plot
% Arnold tongue for each driver

dv1 = unique(dvs(:,1));
dv2 = unique(dvs(:,2));
m1 = length(dv1);
m2 = length(dv2);
indA = 1:m1;
indB = 1:m2;
sync = zeros(m1,m2);

x = linspace(min(dv1),max(dv1),m1);
y = linspace(min(dv2),max(dv2),m2);

figure('name','ArnTongue')
    
for j=1:n
    ind1 = indA(dv1 == dvs(j,1));
    ind2 = indB(dv2 == dvs(j,2));
    sync(ind1,ind2) = S(j,targets(1),targets(2));
end

% plot arnold tongue
imagesc(x,y,sync',[0,1])
set(gca,'YDir','normal')
colormap(jet)
shading(gca,'interp')
colorbar()
title(strcat('Synchrony of Driver ',num2str(targets(1)), ' and Node ',num2str(targets(2))))
xlabel(dependent_vars{1})
ylabel(dependent_vars{2})
%savefig('ArnTongue')

end

