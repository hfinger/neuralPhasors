function [ kernel ] = genSpatialKernel( L, order, sigma, aspectratio, integral, orientation )
%GENSPATIALKERNELS Summary of this function goes here
%   order=-1 generates difference of gaussian (GoD) kernel
%   order=0  generates gaussian kernel
%   order=-1 generates difference of offset gaussian (DooG) kernel
%   order=-1 generates difference of offset gaussian (DooG) kernel


deltax = sigma/2;
kernel=zeros(L);
for xid=1:L
  for yid=1:L
    if mod(L,2)==0
      x = xid-0.5 - L/2;
      y = yid-0.5 - L/2;
    else
      x = xid-1 - L/2;
      y = yid-1 - L/2;
    end
    
    x2=x*cos(orientation) + y*sin(orientation);
    y2=-x*sin(orientation) + y*cos(orientation);
    
    if order==-1
      kernel(xid,yid) = (1+integral)*gauss_1d(x2,sigma)*gauss_1d(y2,sigma*aspectratio) - ...
        gauss_1d(x2,2*sigma)*gauss_1d(y2,2*sigma*aspectratio);
    elseif order==0
      kernel(xid,yid) = gauss_1d(x2,sigma)*gauss_1d(y2,sigma*aspectratio);
    elseif order==1
      kernel(xid,yid) = doog_1d_order1(x2,sigma,deltax,integral)*gauss_1d(y2,sigma*aspectratio);
    elseif order==2
      kernel(xid,yid) = doog_1d_order2(x2,sigma,deltax,integral)*gauss_1d(y2,sigma*aspectratio);
    else
      error('order only up to 2');
    end
  end
end

% if integral~=0
%   %correct the integral in case too small L
%   kernel = integral*kernel/sum(kernel(:));
% end

end

function [out] = gauss_1d(x,sigma)
out = exp( - x .^2 / (2 * sigma^2 ) ) / ( sigma * sqrt(2*pi) );
end

function [out] = doog_1d_order1(x,sigma,deltax,integral)
out = (1+integral) * gauss_1d(x+deltax,sigma) - gauss_1d(x-deltax,sigma);
end

function [out] = doog_1d_order2(x,sigma,deltax,integral)
out = - gauss_1d(x+2*deltax,sigma) + (2+integral) * gauss_1d(x,sigma) - gauss_1d(x-2*deltax,sigma);
end
