function [ drivPosCohsOrdered ] = pltCohOverPhaseOffset( struct, dv, clim, logCoh )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fnames = fieldnames(struct);
dataTmp = struct.(fnames{1});

drivPosCohs = zeros(length(fnames),length(dataTmp.drivPosCoh));
dvs = zeros(length(fnames),1);
for f=1:length(fnames)
    dataTmp = struct.(fnames{f});
    drivPosCohs(f,:) = cell2mat(dataTmp.drivPosCoh);
    dvs(f) = dataTmp.(dv);
end

dvsUnique = unique(dvs);
n = length(dvsUnique);
drivPosCohsOrdered = zeros(n,size(drivPosCohs,2));
yticklabels = cell(1,n);
for i=1:n
    drivPosCohsOrdered(i,:) = mean(drivPosCohs(dvs == dvsUnique(i),:),1);
    yticklabels{i} = num2str(dvsUnique(i));
end

xticklabels = {num2str(0*pi),num2str(0.25*pi),num2str(0.5*pi),num2str(0.75*pi),num2str(1*pi),num2str(1.25*pi),num2str(1.5*pi),num2str(1.75*pi)};
yticks = 1:n;
if logCoh
    imagesc(log(drivPosCohsOrdered))
    title('Log Coherence between driven nodes')
    set(gca,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels)
else
    imagesc(drivPosCohsOrdered)
    title('Coherence between driven nodes')
    set(gca,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels,'CLim',clim)
end
xlabel('Phase Offset in Radians')
ylabel(dv)
colorbar()

end

