function plotKernelDensity( phase, plotdir, densityEstSigma )
%PLOTKERNELDENSITY Summary of this function goes here
%   Detailed explanation goes here


%% Plot kernel densisty estimate:
hfig=figure('Visible','Off');
set(hfig,'MenuBar','none');
clf;
set(hfig,'units','centimeters',...
    'NumberTitle','off','Name','KernelDensity');
pos = get(hfig,'position');
set(hfig,'position',[pos(1:2),3.8,3]);

phi0=linspace(0,2*pi,200);
density=zeros(size(phi0));
for j=1:length(phi0)
    density(j) = densityEstimate( phi0(j) , phase, densityEstSigma );
end
plot(phi0,density);
%xlabel('phase')
%ylabel('kernel density est.')
xlim([0 2*pi]);
set(gca,'XTick',[0 2*pi])
set(gca,'XTickLabel',{'$0$','$2\pi$'})

%save figure:
matlabfrag([plotdir '/KernelDensity']);
end

