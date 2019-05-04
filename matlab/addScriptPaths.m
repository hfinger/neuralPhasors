function addScriptPaths( )
%ADDPATHS Summary of this function goes here
%   Detailed explanation goes here

pathstr = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(pathstr, 'include')))

addpath(fullfile(pathstr));
addpath(fullfile(pathstr, 'classes'));
addpath(fullfile(pathstr, 'methods'));
addpath(fullfile(pathstr, 'methods', 'grid'));
addpath(fullfile(pathstr, 'methods', 'odeSolver'));
addpath(fullfile(pathstr, 'methods', 'plotFeatures'));
addpath(fullfile(pathstr, 'methods', 'connectome'));
addpath(fullfile(pathstr, 'methods', 'connectome_plot'));
addpath(fullfile(pathstr, 'methods', 'metrics'));
addpath(fullfile(pathstr, 'methods', 'metrics','reviewed'));
addpath(fullfile(pathstr, 'ProbtrackClustering'));
addpath(fullfile(pathstr, 'methods', 'connectome_clustering'));



end

