function [variableParams, paramComb] = getVariableParams(paramStruct,combParallel)
% Compute all combinations of variable Parameters:
if isfield(paramStruct,'combParallel') && paramStruct.combParallel
  combParallel=true;
end
[variableParams,paramComb] = getVarSubparams(cell(0),[],paramStruct,cell(0),combParallel);
% paramComb = cell(0);
% for i=1:size(paramCombTmp,2)
%   paramComb{i} = paramCombTmp(:,i);
% end

end

function [variableParams,paramComb] = getVarSubparams(variableParams,paramComb,paramStruct,curVariableParamsDesc,combParallel)
subparamfields = fieldnames(paramStruct);
for i=1:length(subparamfields)
  if length(paramStruct)>1
    
  elseif iscell(paramStruct.(subparamfields{i}))
    %% add this variable param:
    variableParams{end+1,1}=cat(2,curVariableParamsDesc,subparamfields{i}); %#ok<AGROW>
    if isempty(paramComb)
      paramComb=paramStruct.(subparamfields{i});
    else
      if combParallel
        thisParams = paramStruct.(subparamfields{i});
        paramComb(end+1,:)=thisParams; %#ok<AGROW>
      else
        paramComb=combveccell(paramComb,paramStruct.(subparamfields{i}));
      end
    end
  elseif isstruct(paramStruct.(subparamfields{i}))
    %% enter recursive because substructure:
    [variableParams,paramComb] = getVarSubparams(variableParams, ...
      paramComb, ...
      paramStruct.(subparamfields{i}), ...
      cat(2,curVariableParamsDesc,subparamfields{i}), ...
      combParallel);
  end
end
end