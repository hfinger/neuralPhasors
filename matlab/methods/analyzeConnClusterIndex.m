function stats = analyzeConnClusterIndex( W )
%ANALYZECONNECTIVITY Summary of this function goes here

fids=sort(unique(W.f0));

%% format of W:
% neuron(x,y,f1) connects to neuron(x+dx,y+dy,f0)
%% format of conns:
% neuron(x,y,fid) connects to neuron(x+dx,y+dy,toFid)

conns=cell(max(fids),1);
for fid=1:max(fids)
  ids=find(W.f0==fid);
  conns{fid}.dx = -W.dx(ids);
  conns{fid}.dy = -W.dy(ids);
  conns{fid}.toFid = W.f1(ids);
  conns{fid}.w = W.w(ids);
end

for fid=1:max(fids)
  ids=find(W.f1==fid);
  conns{fid}.dx = cat(1,conns{fid}.dx,W.dx(ids));
  conns{fid}.dy = cat(1,conns{fid}.dy,W.dy(ids));
  conns{fid}.toFid = cat(1,conns{fid}.toFid,W.f0(ids));
  conns{fid}.w = cat(1,conns{fid}.w,W.w(ids));
end

stats.clustering_coeffs = zeros(max(fids),1);
stats.meanConnLength = zeros(max(fids),1);
stats.meanConnLengthExc = zeros(max(fids),1);
stats.meanConnLengthInh = zeros(max(fids),1);
for fid=1:max(fids)
  disp(num2str(fid));
  allNgh=conns{fid};
  numConns=0;
  for i=1:length(allNgh.toFid)
    dx=allNgh.dx(i);
    dy=allNgh.dy(i);
    for j=(i+1):length(allNgh.toFid)
      totaldx=dx-allNgh.dx(j);
      totaldy=dy-allNgh.dy(j);
      numConns = numConns + any(conns{allNgh.toFid(i)}.dx==totaldx & conns{allNgh.toFid(i)}.dy==totaldy & conns{allNgh.toFid(i)}.toFid==allNgh.toFid(j));
    end
  end
  stats.clustering_coeffs(fid) = numConns / (length(allNgh.toFid)*(length(allNgh.toFid)-1));
  stats.meanConnLength(fid) = mean(sqrt(conns{fid}.dx.^2 + conns{fid}.dy.^2));
end
stats.meanConnLengthAll = mean(sqrt(W.dx.^2 + W.dy.^2));

end