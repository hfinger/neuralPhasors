function [ F ] = getEmpFC( resort, fTarget, singleHemisphere )
% GETCONNECTOME Function that load connectivity and distance matrix of the
% connectome distinguished into 66 different areas
%
%   Input Parameters:
%
%       resort - if true, areas get resorted according to Marlene's indices
%       p - Order of the norm used to normalize the rows of the
%           connectivity matrix. Set to 0 if no normalization is to be done
%       homotopeScaling - Scalar with which homotopic connections will be
%                         scaled that are added to the connectivity matrix
%       singleHemisphere - if true, only return left hemisphere values
%
%   Output:
%
%       C - 66 x 66 connectivity matrix
%       D - 66 x 66 distance matrix

%% load C and D, resort them and normalize C

% load raw connectivity matrices
path_FCmat = '/media/hofinger/OS/Users/hofinger/phasesim/dti2eeg/eeg_lcmv/COH.mat';

FC = load(path_FCmat);
F = squeeze(mean(mean(mean(FC.FC,1),2),3));
F = F(:,:,fTarget);

% set
F(isnan(F))=0;
F = F + F';

% resort C and D
if resort
    path_ResortIDs = '/media/hofinger/OS/Users/hofinger/phasesim/dti2eeg/sources_plot_order.mat';
    %path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
    resortIds = load(path_ResortIDs);
    resortIds = resortIds.sort_ids;
    F = F(resortIds,:);
    F = F(:,resortIds);
end

if singleHemisphere
    F = F(1:33,1:33);
end

end

