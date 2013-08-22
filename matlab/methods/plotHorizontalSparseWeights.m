function plotHorizontalSparseWeights( inWeightFile )
%PLOTCOVANDCONN Summary of this function goes here
%   Detailed explanation goes here

load(inWeightFile)

%% reorder fids:
fids = W.f0;
orientId=mod(fids,8); orientId(orientId==0)=8;
colorId=ceil(fids/8);
W.f0 = colorId+(orientId-1)*6;

fids = W.f1;
orientId=mod(fids,8); orientId(orientId==0)=8;
colorId=ceil(fids/8);
W.f1 = colorId+(orientId-1)*6;



%%


colormap hsv;
syncIds=find(W.w==-1 );
scatter3(W.dx(syncIds),W.dy(syncIds),W.f1(syncIds),5,W.f0(syncIds))
% scatter3(W.dx(syncIds),W.dy(syncIds),W.f1(syncIds),5,mod(W.f0(syncIds)-24,48))

axis equal
maxr = max(max(abs(W.dx),abs(W.dy)));
xlim([-maxr maxr])
ylim([-maxr maxr])
zlim([0 max(W.f0)])

end