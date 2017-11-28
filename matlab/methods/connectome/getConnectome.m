function [ C,D ] = getConnectome( resort, p, homotopeScaling, singleHemisphere )
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
path_SCmat = '/net/store/nbp/projects/phasesim/databases/avg_SC.mat';
load(path_SCmat);

% resort C and D
if resort
    path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
    load(path_ResortIDs);
    resortIds = [resortIdsMarlene, resortIdsMarlene + 33];
    C = C(resortIds,:);
    C = C(:,resortIds);
    D = D(resortIds,:);
    D = D(:,resortIds);
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
end

end

