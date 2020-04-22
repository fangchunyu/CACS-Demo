function trimmed = trim3d(block, margin) 
%% trim the 3-D block, removing its aliasing on the edge
% Params:
%   margin - [margin_x, margin_y, margin_z] number of pixels to be removed

% % assert length(size(block)) == 3
% % assert length(size(margin))    == 3

trimmed = block(margin(1)+1 : end - margin(1), margin(2)+1 : end - margin(2), margin(3)+1 : end - margin(3));
