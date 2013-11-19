%%  source.avg.filter{grid.inside(indices(i))} has size [3 63]


for i= 1:length(indices) % 66 sources
  for ch=1:63 % 63 kanäle
    svd_s=svd(source.avg.filter{grid.inside(indices(i))}(:,ch)*source.avg.filter{grid.inside(indices(i))}(:,ch)');
    filter_(ch,i)=svd_s(1);
  end
end
source=data;
source.label=sources;
%channel order checked: is the same
for i = 1:length(data.trial)
  source.trial{i} =(data.trial{i}(:,:)'*filter_)';
end