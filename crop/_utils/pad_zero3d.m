function padded = pad_zero3d(stack, padding, mode)
% % wrap a 3-D stack with zeros, 
%   i.e.,padding zeros before and after the content along each dimesion)
%   Params:
%     -padding : a scalar or a 3-elememt vector (3 elements designate numbers of padded zeros along each dimension )
%     -mode : 'wrap' -- padding zeros at both the begining and the end of
%               each dimension;
%             'end' -- padding zeros only at the end of each dimension
        
padding = floor(padding) + 1;
[h, w, d] = size(stack);

if isscalar(padding)
    padding_x = padding, padding_y = padding, padding_z = padding;
else
    if length(padding) == 3
        padding_x = padding(2);
        padding_y = padding(1);
        padding_z = padding(3);
    else error('incorrect format of padding ');
    end
end

if strcmp(mode, 'wrap')
    padded = zeros(h + padding_y*2, w + padding_x*2, d + padding_z*2);
    padded(padding_y : padding_y + h -1, padding_x : padding_x + w - 1, padding_z : padding_z + d - 1) = stack;
elseif strcmp(mode, 'end')
    padded = zeros(h + padding_y, w + padding_x, d + padding_z);
    padded(1 : h, 1 : w, 1 :  d) = stack;
else 
    error('unkonwn mode : %s', mode)
end