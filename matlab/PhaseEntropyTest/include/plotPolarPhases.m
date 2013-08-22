function plotPolarPhases( phase, plotdir )
%PLOTPOLARPHASES Summary of this function goes here
%   Detailed explanation goes here

n=length(phase);

%% Plot phase variables in polar plot:
hfig=figure(51);
clf;
tmp=polar([phase phase]',[zeros(n,1) ones(n,1)]');
set(tmp,'Color','b')
%hide grid:
set(0,'ShowHiddenHandles','on')
set(findobj(findobj(get(gca,'Children'),'Type','Line'),'LineStyle',':'),{'Visible'},{'off'});
%change labels:
dhandles = findall(gcf, 'Type', 'text');
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
    'NumberTitle','off','Name','PolarPhases');
pos = get(hfig,'position');
set(hfig,'position',[pos(1:2),3,3]);
matlabfrag([plotdir '/PolarPhases']);
end

