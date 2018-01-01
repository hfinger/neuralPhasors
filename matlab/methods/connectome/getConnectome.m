function [ C,D,F ] = getConnectome( resort, p, homotopeScaling, singleHemisphere )
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
path_SCmat = '/media/hofinger/OS/Users/hofinger/phasesim/dti2eeg/dti/SC.mat';
path_Distmat = '/media/hofinger/OS/Users/hofinger/phasesim/dti2eeg/dti/dist.mat';
path_FCmat = '/media/hofinger/OS/Users/hofinger/phasesim/dti2eeg/eeg_lcmv/COH.mat';

SC = load(path_SCmat);
Dist = load(path_Distmat);
FC = load(path_FCmat);
C = mean(cat(3,SC.SC{:}),3);
D = mean(cat(3,Dist.dist{:}),3);
F = squeeze(mean(mean(mean(FC.FC,1),2),3));
F = F(:,:,9);

% set
C(isnan(C))=0;
D(isnan(D))=0;
F(isnan(F))=0;
C = C + C';
D = D + D';
F = F + F';

% resort C and D
if resort
    path_ResortIDs = '/media/hofinger/OS/Users/hofinger/phasesim/dti2eeg/sources_plot_order.mat';
    %path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
    resortIds = load(path_ResortIDs);
    resortIds = resortIds.sort_ids;
    C = C(resortIds,:);
    C = C(:,resortIds);
    D = D(resortIds,:);
    D = D(:,resortIds);
    F = F(resortIds,:);
    F = F(:,resortIds);
end

% add homotopic connections to C
C = C + homotopeScaling * diag(ones(size(C,1)/2,1),size(C,1)/2) + homotopeScaling * diag(ones(size(C,1)/2,1),-size(C,1)/2);

% normalize C
if p > 0
    C = bsxfun(@rdivide,C,sum(C.^p,2).^(1/p));
end

if singleHemisphere
    C = C(1:33,1:33);
    D = D(1:33,1:33);
    F = F(1:33,1:33);
end

end

