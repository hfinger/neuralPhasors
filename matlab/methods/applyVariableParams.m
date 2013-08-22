function [paramStruct] = applyVariableParams(paramStruct,variableParams,paramComb,jobid)

for i=1:length(variableParams)
  paramStruct = applyOneVarSubparams(paramStruct,variableParams{i},paramComb{i,jobid});
end

end

function [paramStruct] = applyOneVarSubparams(paramStruct,variableParams,setValue)

fname = variableParams{1};
variableParams(1)=[];
if isempty(variableParams)
  paramStruct.(fname) = setValue;
else
  paramStruct.(fname) = applyOneVarSubparams(paramStruct.(fname),variableParams,setValue);
end
end