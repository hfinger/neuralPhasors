function [ L, H2, C2, H2mean, C2mean, Lmean ] = colorStatsOfPatches( W )
%COLORSTATSOFcolorsES Summary of this function goes here
%   W should be 4 dimensional with dimensions [dx,dy,color,patchid]

colors = permute(W, [1 2 4 3]);
colors = reshape(colors,[size(colors,1)*size(colors,2)*size(colors,3) size(colors,4)]);

  L = mean(colors,2); %luminance
  colorsEqualLum = bsxfun( @minus, colors, L);
  S = sqrt(sum(colorsEqualLum.^2,2)); %saturation
  
  
  %from wikipedia: http://en.wikipedia.org/wiki/HSL_and_HSV
  alpha = 0.5*(2*colors(:,1) - colors(:,2) - colors(:,3));
  beta = 0.5*sqrt(3)*(colors(:,2) - colors(:,3));
  H2 = atan2(beta,alpha);
  C2 = sqrt(alpha.^2 + beta.^2);
  
  
  
  H2 = reshape(H2,[size(W,1) size(W,2) size(W,4)]);
  C2 = reshape(C2,[size(W,1) size(W,2) size(W,4)]);
  L = reshape(L,[size(W,1) size(W,2) size(W,4)]);
  
  H2mean = sum(sum(H2,1),2)/(size(W,1)*size(W,2));
  C2mean = sum(sum(C2,1),2)/(size(W,1)*size(W,2));
  Lmean = sum(sum(L,1),2)/(size(W,1)*size(W,2));
  
  H2mean = mod(H2mean(:),pi);
  C2mean = C2mean(:);
  Lmean = Lmean(:);

end

