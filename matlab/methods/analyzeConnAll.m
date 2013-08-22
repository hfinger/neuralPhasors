function analyzeConnAll( weightFolder, weightFile )
%ANALYZECONNALL Summary of this function goes here
%   Detailed explanation goes here

disp('start all conn')
W = load(fullfile(weightFolder,weightFile));
W = W.W;

[ stats ] = analyzeConnSmallWorldIndex( W );
save(fullfile(weightFolder,'stats.mat'),'stats','W')

disp('start only sync conn')
W = load(fullfile(weightFolder,weightFile));
W = W.W;
syncIds = find(W.w>0);
W.dx = W.dx(syncIds);
W.dy = W.dy(syncIds);
W.f0 = W.f0(syncIds);
W.f1 = W.f1(syncIds);
W.w = W.w(syncIds);
[ stats ] = analyzeConnSmallWorldIndex( W );
save(fullfile(weightFolder,'statsOnlySyncConn.mat'),'stats','W')

disp('start only desync conn')
W = load(fullfile(weightFolder,weightFile));
W = W.W;
desyncIds = find(W.w<0);
W.dx = W.dx(desyncIds);
W.dy = W.dy(desyncIds);
W.f0 = W.f0(desyncIds);
W.f1 = W.f1(desyncIds);
W.w = W.w(desyncIds);
[ stats ] = analyzeConnSmallWorldIndex( W );
save(fullfile(weightFolder,'statsOnlyDesyncConn.mat'),'stats','W')

end

