function [ matOut ] = normMat( mat )
%NORMMAT Normalizes a matrix in both directions using regression
%   Detailed explanation goes here

% Run Test:
% dim1vars = rand(4,1);
% dim2vars = rand(4,1);
% matIn = bsxfun(@plus, dim1vars, dim2vars');
% matIn = matIn + 0.3*eye(size(matIn));
% [ matOut ] = normMat( matIn );
% figure(1); imagesc( matIn )
% figure(2); imagesc( matOut )

dim1Group = repmat(1:size(mat,1),[size(mat,2) 1])';
dim2Group = repmat(1:size(mat,2),[size(mat,1) 1]);

predictors = dummyvar([dim1Group(:) dim2Group(:)]);
[~,~,r] = regress(mat(:),predictors(:,1:end-1));
matOut = reshape(r,size(mat));


end

