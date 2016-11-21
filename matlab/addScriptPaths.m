function addScriptPaths( )
%ADDPATHS Summary of this function goes here
%   Detailed explanation goes here

pathstr = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(pathstr, 'include')))

addpath(fullfile(pathstr));
addpath(fullfile(pathstr, 'classes'));
% addpath(fullfile(pathstr, 'include'));
% addpath(fullfile(pathstr, 'include', 'append_pdfs'));
% addpath(fullfile(pathstr, 'include', 'CircStat2012a'));
% addpath(fullfile(pathstr, 'include', 'combveccell'));
% addpath(fullfile(pathstr, 'include', 'discretesample'));
% addpath(fullfile(pathstr, 'include', 'BCT'));
% addpath(fullfile(pathstr, 'include', 'barweb'));
% addpath(fullfile(pathstr, 'include', 'export_fig'));
% addpath(fullfile(pathstr, 'include', 'FastICA_25'));
% addpath(fullfile(pathstr, 'include', 'fdr_bh'));
% addpath(fullfile(pathstr, 'include', 'fprintf_r'));
% addpath(fullfile(pathstr, 'include', 'CircStat2012a'));
% addpath(fullfile(pathstr, 'include', 'cline'));
% addpath(fullfile(pathstr, 'include', 'gifti-1.4'));
% addpath(fullfile(pathstr, 'include', 'grs2rgb'));
% addpath(fullfile(pathstr, 'include', 'matlabfrag'));
% addpath(fullfile(pathstr, 'include', 'matlab2tikz_0.3.2','helpers'));
% addpath(fullfile(pathstr, 'include', 'matlab2tikz_0.3.2','src'));
% addpath(fullfile(pathstr, 'include', 'minFunc_2012', 'minFunc'));
% addpath(fullfile(pathstr, 'include', 'minFunc_2012', 'minFunc', 'compiled'));
% addpath(fullfile(pathstr, 'include', 'plot2svg_20120915'));
% addpath(fullfile(pathstr, 'include', 'safeLoad'));
% addpath(fullfile(pathstr, 'include', 'shadedErrorBar'));
% addpath(fullfile(pathstr, 'include', 'sigstar'));
% addpath(fullfile(pathstr, 'include', 'sfigure'));
% addpath(fullfile(pathstr, 'include', 'panel_2.10'));
% addpath(fullfile(pathstr, 'include', 'subaxis'));
% addpath(fullfile(pathstr, 'include', 'TLS'));
addpath(fullfile(pathstr, 'methods'));
addpath(fullfile(pathstr, 'methods', 'grid'));
addpath(fullfile(pathstr, 'methods', 'odeSolver'));
addpath(fullfile(pathstr, 'methods', 'plotFeatures'));
addpath(fullfile(pathstr, 'methods', 'connectome'));
addpath(fullfile(pathstr, 'methods', 'connectome_plot'));
addpath(fullfile(pathstr, 'methods', 'metrics'));
addpath(fullfile(pathstr, 'methods', 'metrics','reviewed'));
addpath(fullfile(pathstr, 'ProbtrackClustering'));

end

