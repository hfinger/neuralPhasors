function [ x,y,z ] = combineClusterSizeAndEntrop( a, b, c, d )
%PLOTCLUSTERSIZEANDENTROP Summary of this function goes here
%   Detailed explanation goes here

tmp=evalin('base','out;');

out=tmp{d};

x1=a*out.data.entropy.^b;
y1=out.data.logClusterSize.^c;
z=(x1+y1)';

x=out.xData;
y=out.yData;
end

