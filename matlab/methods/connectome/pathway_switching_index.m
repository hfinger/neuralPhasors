function [switching_index] = pathway_switching_index(PAs_of_collected_paths1,PAs_of_collected_paths2)
%PATHWAY_SWITCHING_INDEX Summary of this function goes here
%   Detailed explanation goes here

PA1_minus_PA2 = squeeze(PAs_of_collected_paths1(1:16, :) - PAs_of_collected_paths2(1:16, :));
PA2_minus_PA1 = squeeze(PAs_of_collected_paths2(1:16, :) - PAs_of_collected_paths1(1:16, :));
max_PA1_minus_PA2 = max(PA1_minus_PA2, [], 1);
max_PA2_minus_PA1 = max(PA2_minus_PA1, [], 1);
switching_index = max_PA1_minus_PA2 .* max_PA2_minus_PA1;

%switching_index = sign(switching_index) .* sqrt(abs(switching_index));
switching_index = switching_index ./ sqrt(abs(switching_index));

end

