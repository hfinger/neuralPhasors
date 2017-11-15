function [ CohsOrdered, uniqueDVs ] = pltCohOver2DVs( struct, dvs, targets, pltCoh, clim, logCoh, meshplt )
%PLTCOHOVER2DVS Function that plots the coherence over two dependent
%variables
%   Input Parameters:
%       
%       struct - structure that includes 1 field for each simulation. Each
%                field needs to have a subfield called FC where the
%                coherence value is stored
%       dvs - cell with 2 strings refering to the dependent variables. Need
%             to be a subfield of each field of struct
%       pltCoh - if false, no plot is created
%       clim - vector with upper and lower bound for the color plot
%       logCoh - if true, log coherence will be plotted
%       meshplt - if true, surface plot will be created. Else imagesc
%
%   Output:
%
%       drivPosCohsOrdered - 2-dimensional array with the coherence values
%                            ordered by the 2 dvs
%       uniqueDVs = 2-dimensional array with the sorted, unique values of
%                   both dvs

%% extract coherence and dvs from struct

% get fieldnames and initialize stuff
fnames = fieldnames(struct);
dataTmp = struct.(fnames{1});
Cohs = zeros(length(fnames),1);
dv_coll = zeros(length(fnames),2);

for f=1:length(fnames)
    
    dataTmp = struct.(fnames{f});
    
    % extract coherence
    Coh = dataTmp.FC;
    if iscell(Coh)
        Cohs(f) = Coh{1,1}(targets(1),targets(2));
    else
        Cohs(f) = Coh(targets(1),targets(2));
    end
    
    % extract dependent variables
    dv_coll(f,1) = dataTmp.(dvs{1,1});
    dv_coll(f,2) = dataTmp.(dvs{1,2});
    
end

%% order the coherences according to dvs and collect plot labels

% get unique values of both dvs
dv1 = unique(dv_coll(:,1));
dv2 = flipud(unique(dv_coll(:,2)));

% initialize plot label vectors
CohsOrdered = zeros(length(dv1),length(dv2));
xtickTargets = 1:ceil(length(dv1)/8):length(dv1);
xticks = zeros(1,length(xtickTargets));
yticks = 1:length(dv2);
xticklabels = cell(1,length(xticks));
yticklabels = cell(1,length(yticks));

% loop over both dvs
n = 1;
for i=1:length(dv1)
    
    % get index for first dv
    idx1 = dv_coll(:,1) == dv1(i);
    
    for j=1:length(dv2)
        
        % get index for second dv
        idx2 = dv_coll(:,2) == dv2(j);
        
        % extract coherence value at indice position
        CohsOrdered(i,j) = mean(Cohs(idx1.*idx2 == 1,:),1);
        
        % store y axis plot labels
        if i == 1
            yticklabels{j} = num2str(dv2(j));
        end
        
    end
    
    % store x axis plot labels
    if sum(xtickTargets == i) > 0
        xticks(n) = dv1(i);
        xticklabels{n} = num2str(dv1(i));
        n = n+1;
    end
    
end

% store dv values
uniqueDVs{1} = dv1;
uniqueDVs{2} = dv2;

%% Plotting

if pltCoh
    if logCoh
        drivPosCohsOrdered = log(CohsOrdered);
    end

    figure()
    if meshplt
        h = meshz(dv1, dv2, CohsOrdered');
        xlim([min(dv1) max(dv1)])
        set(gca,'XTick',xticks)
        set(gca,'XTickLabel',xticklabels)
        zlabel('Coherence')
        ylim([min(dv2) max(dv2)])
        set(gca,'fontsize',14)
        set(gca,'linewidth',2)
        set(h, 'faceColor', 'flat')
        set(h, 'edgeColor', 'k')
        axis vis3d
        az = 55;
        el = 40;
        view(az, el);
    else
        imagesc(CohsOrdered')
        set(gca,'XTick',xtickTargets,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels,'CLim',clim)
        colorbar()
    end
    title('Coherence between driven nodes')
    xlabel(dvs{1,1})
    ylabel(dvs{1,2})
end

end

