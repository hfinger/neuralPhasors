function [ dvs_sorted ] = plotConnMats( dataStruct, dependent_vars, coherence_measure )
%PLOTCONNMATS Function that creates subplots of coherence measure matrixes
% ordered by the unique values the parameters in dependent_vars can take on 
%   Input Parameters:
%       dataStruct - structure with 1 field for each file
%       dependent_vars - 1 or 2 dependent variables (cell of strings)
%       coherence_measure - SE for Shannon Entropy or Coherence (string)
%   Returns:
%       dvs_sorted - array including the values of the dvs at each subplot

fnames = fieldnames(dataStruct);
n = length(fnames);
if ~iscell(dependent_vars)
    dvs_tmp = dependent_vars;
    dependent_vars = cell({});
    dependent_vars{1} = dvs_tmp;
end
m = length(dependent_vars);
if m > 2
    error('Function can handle 2 dependent variables at maximum')
end
dvs = zeros(n,m);

for f = 1:n
    
    data = dataStruct.(fnames{f});
    
    FC{f} = data.(coherence_measure);
    
    for d = 1:m
        dvs_tmp = data.(dependent_vars{d});
        if length(dvs_tmp) > 1
            dvs_tmp = dvs_tmp(1);
        end
        dvs(f,d) = dvs_tmp; 
    end
    
end

for i = 1:m  
    dvs_sorted{i} = sort(unique(dvs(:,i)));
end

figure('name','ConnMats','units','normalized','position',[.1,.1,.7,.7],'outerposition',[0,0,1,1])
idx = 1;
for k = 1:length(dvs_sorted{1})
    for l = 1:length(dvs_sorted{2})
        subplot(length(dvs_sorted{1}),length(dvs_sorted{2}),idx)
        ind1 = dvs(:,1) == dvs_sorted{1}(k);
        ind2 = dvs(:,2) == dvs_sorted{2}(l);
        ind = ind1 == 1 & ind2 == 1;
        imshow(FC{ind},[0 1],'InitialMagnification', 'fit')
        colormap('jet')
        title(strcat(dependent_vars{1},' = ',num2str(dvs_sorted{1}(k)),' , ',dependent_vars{2},' = ',num2str(dvs_sorted{2}(l))))
        if l == length(dvs_sorted{2})
            colorbar()
        end
        idx = idx + 1;
    end
end
savefig('ConnMats')
end

