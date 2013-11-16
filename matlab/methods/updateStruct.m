function [ oldStruct ] = updateStruct( oldStruct,newStruct )
%update a structure with fields of newStruct:
fnames=fieldnames(newStruct);

if length(newStruct)>1
  
  %special case of structure array:
  
%   for j=1:length(oldStruct)
%     for i=1:length(fnames)
%         if isstruct(newStruct(j).(fnames{i})) && isfield(oldStruct(j),fnames{i}) && isstruct(oldStruct(j).(fnames{i}))
%           oldStruct(j).(fnames{i}) = updateStruct(oldStruct(j).(fnames{i}),newStruct(j).(fnames{i}));
%         else
%           oldStruct(j).(fnames{i}) = newStruct(j).(fnames{i});
%         end
%     end
%   end
  oldStruct = newStruct;
  disp('WARNING: structure array was not updated with subfield but just assigned... in updateStruct(oldStruct,newStruct)')
  
else
  
  %normal case of simple structure:
  
  for i=1:length(fnames)
      if isstruct(newStruct.(fnames{i})) && isfield(oldStruct,fnames{i}) && isstruct(oldStruct.(fnames{i}))
        oldStruct.(fnames{i}) = updateStruct(oldStruct.(fnames{i}),newStruct.(fnames{i}));
      else
        oldStruct.(fnames{i}) = newStruct.(fnames{i});
      end
  end
  
end

end

