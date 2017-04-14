function [ dataStruct ] = getArgsFromFiles( resultPath, identifier, coherence_measure, incData, coh_index, varargin )
%GETARGSFROMFILES Reads in all data from resultpath including identifier in
%their name, extracts the coherence_measure and simulation parameters
%defined by varargin and puts them into a struct
%   Input Parameters:
%       resultPath - Full directory to files (string)
%       identifier - characters included in all target files, e.g. *.mat (string)
%       coherence_measure - SE for Shannon Entropy or Coherence (string)
%       varargin - variable number of desired simulation parameters (string)
%   Returns:
%       dataStruct - struct with 1 field for each file

files = dir(fullfile(resultPath,identifier));
fnames = {files.name};
n = length(fnames);

for f=1:n
    data = load(strcat(resultPath, '/', fnames{f}));
    struct_tmp.Coherence = data.(coherence_measure);
    if coh_index > 0
        struct_tmp.Coherence = struct_tmp.Coherence(:,:,coh_index);
    end
    if incData
        struct_tmp.Y = data.simResult.Y;
    end
    for v=1:length(varargin)
        struct_tmp.(varargin{v}) = data.simResult.sim.(varargin{v});
    end
    dataStruct.(fnames{f}(fliplr(end-4:-1:1))) = struct_tmp;
end

end
