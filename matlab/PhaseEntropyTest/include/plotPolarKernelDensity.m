function plotPolarKernelDensity( phase, plotdir, densityEstSigma )
%PLOTKERNELDENSITY Summary of this function goes here
%   Detailed explanation goes here


%% Plot kernel densisty estimate:
hfig=figure('Visible','Off');
set(hfig,'MenuBar','none');
clf;
phi0=linspace(0,2*pi,200);
density=zeros(size(phi0));
for j=1:length(phi0)
    density(j) = densityEstimate( phi0(j) , phase, densityEstSigma );
end
tmp=polar(phi0,density);
set(tmp,'Color','b')
%hide grid:
%set(0,'ShowHiddenHandles','on')
%set(findobj(findobj(get(gca,'Children'),'Type','Line'),'LineStyle',':'),{'Visible'},{'off'});
%change labels:
dhandles = findall(gcf, 'Type', 'text');%, 'HorizontalAlignment', 'center');
for K = 1:length(dhandles)
    deg = str2double( get(dhandles(K), 'String') );
    if deg==0
        set(dhandles(K),'String','$0$');
    elseif deg==90
        set(dhandles(K),'String','$\pi/2$');
    elseif deg==180
        set(dhandles(K),'String','$\pi$');
    elseif deg==270
        set(dhandles(K),'String','$3\pi/2$');
    else
        set(dhandles(K),'String','');
    end
end

%save figure:
set(hfig,'units','centimeters',...
    'NumberTitle','off','Name','PolarKernelDensity');
pos = get(hfig,'position');
set(hfig,'position',[pos(1:2),3.8,3.8]);
matlabfrag([plotdir '/PolarKernelDensity']);
end

