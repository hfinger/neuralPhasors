function [ newStructs ] = newStruct( data, dv, dv_extract, extract_idx )
%NEWSTRUCT Function that creates new structures, one for each value of the
%given variable
% 
%   Input Parameters:
%       data - old structure
%       dv - variable, which has to be a field of each entry in struct
% 
%   Output:
%       newStructs - structure with entries being one structure for each
%       value of dv

%%

fnames = fieldnames(data);
newStructs = struct();
dvs = cell(0);

if ~isempty(dv_extract)
    dv_extract_tmp = data.(fnames{1}).(dv_extract);
    if length(size(dv_extract_tmp)) == 2
        extract_idx = sub2ind(size(dv_extract_tmp),extract_idx(1),extract_idx(2));
    end
end
        
for f=1:length(fnames)
    
    data_tmp = data.(fnames{f});
    dv_tmp = data_tmp.(dv);
    
    if isempty(dvs)
        i = 1;
        dvs{i} = dv_tmp;
    else
        check = false;
        for j=1:length(dvs)
            if dvs{j} == dv_tmp
                i = j;
                check = true;
                break
            end
        end
        if ~check
            i = length(dvs) + 1;
            dvs{i} = dv_tmp;
        end
    end
    
    if ~isempty(dv_extract)
        dv_extract_tmp = data_tmp.(dv_extract);
        data_tmp.(dv_extract) = dv_extract_tmp(extract_idx);
    end
    
    newStructs.(['newStruct_' num2str(i)]).(fnames{f}) = data_tmp;
    
end

end

