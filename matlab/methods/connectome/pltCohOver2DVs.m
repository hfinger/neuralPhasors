function [ drivPosCohsOrdered, uniqueDVs ] = pltCohOver2DVs( struct, dvs, pltCoh, clim, logCoh, meshplt )
%PLTCOHOVER2DVS Function that plots the coherence over two dependent
%variables
%   Input Parameters:
%       
%       struct - structure that includes 1 field for each simulation. Each
%                field needs to have a subfield called drivPosCoh where the
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
fnames = fieldnames(struct);
dataTmp = struct.(fnames{1});

drivPosCohs = zeros(length(fnames),length(dataTmp.drivPosCoh));
dv_coll = zeros(length(fnames),2);
for f=1:length(fnames)
    dataTmp = struct.(fnames{f});
    drivPosCohs(f,:) = cell2mat(dataTmp.drivPosCoh);
    dv_coll(f,1) = dataTmp.(dvs{1,1});
    dv_coll(f,2) = dataTmp.(dvs{1,2});
end

%% order the coherences according to dvs and collect plot labels

% get unique values of both dvs
dv1 = unique(dv_coll(:,1));
dv2 = unique(dv_coll(:,2));

% initialize plot label vectors
drivPosCohsOrdered = zeros(length(dv1),length(dv2));
xtickTargets = 1:ceil(length(dv1)/8):length(dv1);
xticks = zeros(1,length(xtickTargets));
yticks = 1:length(dv2);
xticklabels = cell(1,length(xticks));
yticklabels = cell(1,length(yticks));

% loop over both dvs
n = 1;
for i=1:length(dv1)
    idx1 = dv_coll(:,1) == dv1(i);
    for j=1:length(dv2)
        idx2 = dv_coll(:,2) == dv2(j);
        drivPosCohsOrdered(i,j) = mean(drivPosCohs(idx1.*idx2 == 1,:),1);
        if i == 1
            yticklabels{j} = num2str(dv2(j));
        end
    end
    if sum(xtickTargets == i) > 0
        xticks(n) = dv1(i);
        xticklabels{n} = num2str(dv1(i));
        n = n+1;
    end
end

if pltCoh
    if logCoh
        drivPosCohsOrdered = log(drivPosCohsOrdered);
    end

    figure()
    if meshplt
        h = meshz(dv1, dv2, drivPosCohsOrdered');
        xlim([0 2*pi])
        set(gca,'XTick',xticks)
        set(gca,'XTickLabel',xticklabels)
        zlabel('Coherence')
        ylim([0 0.027])
        set(gca,'fontsize',14)
        set(gca,'linewidth',2)
        set(h, 'faceColor', 'flat')
        set(h, 'edgeColor', 'k')
        axis vis3d
        az = 55;
        el = 40;
        view(az, el);
    else
        imagesc(drivPosCohsOrdered')
        set(gca,'XTick',xticks,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels,'CLim',clim)
        colorbar()
    end
    title('Coherence between driven nodes')
    xlabel(dvs{1,1})
    ylabel(dvs{1,2})
end

uniqueDVs{1} = dv1;
uniqueDVs{2} = dv2;

end

