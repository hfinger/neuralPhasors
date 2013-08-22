function plotTestcaseSurf( out, savedir )
%PLOTTESTCASEOPTSURF Summary of this function goes here
%   Detailed explanation goes here
addpath('include')

fnames=fieldnames(out.data);
[Selection,ok] = listdlg('PromptString','Select a variable:',...
  'SelectionMode','multiple',...
  'ListString',fnames);
for i=1:length(Selection)
  hfig=figure;
  set(hfig,'units','centimeters',...
    'NumberTitle','off','Name',fnames{Selection(i)});
  pos = get(hfig,'position');
  set(hfig,'position',[pos(1:2),13,9]);
  
  surf(out.xData,out.yData,real(out.data.(fnames{Selection(i)}))')
  ylabel(out.yLabel)
  xlabel(out.xLabel)
  zlabel(fnames{Selection(i)})
  xlim([min(out.xData),max(out.xData)])
  ylim([min(out.yData),max(out.yData)])
  
  colormap winter;
  if nargin>1
    set(gcf, 'color', 'w');
    matlabfrag([savedir '/' fnames{Selection(i)}],'renderer','opengl');
    set(hfig,'units','pixels');
    pos = get(hfig,'position');
    set(hfig,'position',[pos(1:2),13*40,9*40]);
    drawnow;
    I = getframe(gcf);
    imwrite(I.cdata, [savedir '/' fnames{Selection(i)} '.png'],'png');
  end
end

end

