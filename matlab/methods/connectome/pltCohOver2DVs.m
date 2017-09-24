function [ drivPosCohsOrdered ] = pltCohOver2DVs( struct, dvs, pltCoh, clim, logCoh, meshplt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

dv1 = unique(dv_coll(:,1));
dv2 = unique(dv_coll(:,2));
drivPosCohsOrdered = zeros(length(dv1),length(dv2));
xticks = 1:ceil(length(dv1)/8):length(dv1);
yticks = 1:length(dv2);
xticklabels = cell(1,length(xticks));
yticklabels = cell(1,length(yticks));
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
    if sum(xticks == i) > 0
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
        %set(gca,'XTick',xticks)
        %set(gca,'XTickLabel',xticklabels)
        xlim([0 2*pi])
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

end

