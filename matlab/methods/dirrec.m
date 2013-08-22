function [ paths, filenames ] = dirrec( parent , fileRegexp )
%DIRREC Summary of this function goes here
%   Detailed explanation goes here

paths = cell(0);
filenames = cell(0);

% list all sub folders and files:
listing = dir(parent);

% list all files in parent folder:
for i=1:length(listing)
  if listing(i).isdir
    if ~strcmp(listing(i).name,'.') && ~strcmp(listing(i).name,'..')
      [ pathsTmp, filenamesTmp ] = dirrec( fullfile(parent,listing(i).name) , fileRegexp );
      paths(end+1:end+length(pathsTmp)) = pathsTmp;
      filenames(end+1:end+length(filenamesTmp)) = filenamesTmp;
    end
  else
    indexMatches = regexp(listing(i).name, fileRegexp);
    if indexMatches
      paths{end+1} = parent;
      filenames{end+1} = listing(i).name;
    end
  end
end

