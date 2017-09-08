function [ hg ] = plotColoredBrain( cData, resort, transparency )
%PLOTCOLOREDBRAIN Function that colors 66 areas of brain according to cData
%   Input Parameters:
%       cData - vector of length 66 that is used to color 66 brain regions
%       resort - If true, rearange the order of cData
%                (reverses order achieved by resortIdsMarlene)
%   Output:
%       hg - figure handle

%% load brain surface and mapping to 66 distinct regions (ROIs)
addpath('/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e')
subject = 3;
%savefolder = strcat('Coherence_Stimulated_subj',num2str(subject));
%mkdir(savefolder);
subjectStr = num2str(subject,'%02u');

g = gifti(['ca' subjectStr '_1_structcortex_8196.surf.gii']);
regionmapping = importdata([subjectStr '_regionmapping.txt']);

%% use regionmapping and cData to get coloring for each ROI
cdataPerROI = [cData; zeros(1)];
regionmappingPlusOne = regionmapping;
regionmappingPlusOne(regionmappingPlusOne==0) = length(cdataPerROI);

% resort cData if necessary
if resort 
    path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
    load(path_ResortIDs)
    resortIDs = [resortIdsMarlene, resortIdsMarlene + 33, length(cdataPerROI)];
    idx = 1:length(cdataPerROI);
    for i=idx
        idx(resortIDs(i)) = i;
    end
    cdataPerROI = cdataPerROI(idx);
end

cdataPerVertex.cdata = cdataPerROI(regionmappingPlusOne,:);

%% plot surface with ROI-colors:
figure();
clf; 
hg=plot(g,cdataPerVertex);
colormap('cool')
colorbar()
shading flat
if transparency
    alpha(cdataPerVertex.cdata);
end
view(-9.8705,72.0327);

end

