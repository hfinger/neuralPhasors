function   SC = normGraph(SC, avg_roi_size,ROInorm, ROWnorm, sparse)
%% make connectivity of SC sparse                                  do this before or after matrix normalizations?

% values of sparse = [0.6, 0.7, 0.75[, half-open interval (no sign. RC for boundaries)
% Sporns et al (2007) deterministic tracking SC has sparseness 0.7 to 0.75
% Messe et al (2014) probabilisitc tracking SC is thresholded at .001, sparseness?
% Goni et al (2014) ...
% Zamora-Lopez et al (2010) cat cortex data with sparseness 0.7

if (sparse >= 0) && (sparse <= 1)
  [~,sIdx] = sort(SC(:),'ascend');
  mIdx = sIdx(1:int32(sparse*numel(SC)));
  SC(mIdx) = 0;                                                           % set elements to 0 in ascending order,
  [d1,d2] = ind2sub(size(SC), mIdx);                                      % get indices of transpose elements
  SC(sub2ind(size(SC),d2,d1)) = 0;                                        % set to 0 as well to conserve symmetry
end

effSparse = sum(sum(SC == 0)) / numel(SC);                                % determine effective sparseness of resulting graph

%% norm connections by ROI sizes to get connectivity density XXX

if strcmp(ROInorm,'ROIprd')
  ROIprd = avg_roi_size*avg_roi_size';                                    % normalize by product of ROI sizes
  SC = SC ./ ROIprd;                                                      % RCsignif more sensitive to sparseness in normSum than in normROI
elseif strcmp(ROInorm,'ROIsum')
  ROIsum = bsxfun(@plus,avg_roi_size,avg_roi_size');                      % normalize by sum of ROI sizes
  SC = SC ./ ROIsum;
end

if ROWnorm
  % norm rows, i.e. input, to sum(CIJ) == 1:                              % remove this from connectivity data
  SC = bsxfun(@rdivide,SC,sum(SC,2))';                                    % and make it part of the model?
end