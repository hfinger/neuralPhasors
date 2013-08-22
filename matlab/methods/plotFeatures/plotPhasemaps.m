function plotPhasemaps( act, phase )
%PLOTPHASEMAPS Summary of this function goes here
%   Detailed explanation goes here

d1=size(act,1);
d2=size(act,2);
d3=size(act,3);


phase = mod(phase,2*pi)-pi;
phaseReshaped = reshape(permute(reshape(phase,[d1 d2 sqrt(d3) sqrt(d3)]),[1 3 2 4]),[d1*sqrt(d3) d2*sqrt(d3)]);
actReshaped = reshape(permute(reshape(act,[d1 d2 sqrt(d3) sqrt(d3)]),[1 3 2 4]),[d1*sqrt(d3) d2*sqrt(d3)]);

figure(1);
imagesc(actReshaped); 
axis image;
colormap hot;
colorbar;
title('activity')

figure(2); 
imagesc(phaseReshaped,[-pi pi]); 
axis image;
colormap hsv;
colorbar;
title('phase')

figure(3); 
h = imagesc(phaseReshaped,[-pi pi]); 
set(h, 'AlphaData', actReshaped);
alim([min(actReshaped(:)) max(actReshaped(:))]);
axis image;
colormap hsv;
colorbar;
title('phase with activity as alpha channel')

figure(4); 
imagesc(angle(sum(act.*exp(1i * phase),3)),[-pi pi]); 
axis image;
colormap hsv;
colorbar;
title('phase weighted by activity')

end

