function [ missingFiles, minFile, maxFile ] = findMissingFiles( path )
%FINDMISSINGFILES Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
  path = pwd;
end

listing = dir(path);
files = {listing.name};
numbers = regexp(files,'\d+','match');
numbers = cellfun(@(x) str2double(x), numbers, 'UniformOutput', false);
numbers = cell2mat(numbers);

minFile = min(numbers);
maxFile = max(numbers);

missingFiles = setdiff(minFile:maxFile,numbers);

end

