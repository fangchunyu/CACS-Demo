function grid_dim = crop3d_fn(new_crop_size, new_overlap, save_all, margin_mode, varargin)
%% crop small 3-D stacks into training data
%  Params:
%   -save_all: boolean, if false, blocks with pixel value sum below
%       varargin{1} or variance below varargin{2} will be discarded.
%
%   -margin_mode, int in [1, 2, 3] that indicates how to deal with the
%        margin pixels of the 3-D stack (when the margin pixels is not enough as a block)
%
%       1: discard margin pixels
%       2: padding zeros to fill a block. usually used in validation data cropping.
%       3: cover_all: shift the cropping window backwards to make more overlap so that the
%           margin will be covered. useful in training data cropping.




path(path, '_utils')

%crop_size = [100, 100, 20];

%overlap = [0.5, 0.5, 0.5];

%pixel_threshold = 4.15e6;
%var_threshold   = 2e4;

% % for 16bit_bessel brain_block of 80*80*64
% pixel_threshold = 1.5e9;
% var_threshold   = 2.54e5;

% % for bessel brain_block of 256*256*64:
% pixel_threshold = 12e7;
% var_threshold   = 1e2;

% % for bessel brain_block of 160*160*80:
% pixel_threshold = 1.5e7;
% var_threshold   = 5e1;

% % for bessel brain_block of 40*40*40:
% pixel_threshold = 2.5e6;
% var_threshold   = 3.5e1;

% % for celegans of 160*160*64
%pixel_threshold = 17e5;
%var_threshold   = 5e1;

% % for celegans of 80*80*64
% pixel_threshold = 4e5;
% var_threshold   = 5e1;

% % for 3T3_488 of 80*80*80
%pixel_threshold = 6e5;
%var_threshold   = 1.5e1;

% % for 3T3_405(16bit) of 80*80*80
%pixel_threshold = 5e8;
%var_threshold   = 6.7e5;

if ~save_all
    pixel_threshold = varargin{1};
    var_threshold   = varargin{2};
end

padding = false;
cover_all = false;

if margin_mode == 2
    padding = true;
    cover_all = false;
elseif margin_mode == 3
    padding = false;
    cover_all = true;
end


block_bias = 0;


% if ~exist('new_overlap', 'var')
%     new_overlap = [0., 0., 0.];
% end

overlap_px = floor(new_overlap .* new_crop_size);

if exist('preferences.mat', 'file')
    load preferences.mat
    filepath = crop_filepath;
else
    filepath = pwd;
end
[file_name,  filepath] = uigetfile({ '*.tif'; '*.*'}, 'MultiSelect', 'on', 'Select stack to crop...', filepath);

if ~iscell(file_name)
    if file_name == 0
        return
    end
    file_name = {file_name};
end


sampleFile = fullfile(filepath, file_name{1});
sampleInfo = imfinfo(sampleFile);
bitDepth = sampleInfo.BitDepth

file_num = length(file_name);
cropped_num = 0;
%save_dir = sprintf('\\cropped%dX%dX%d_overlap%.1f-%.1f-%.1f', new_crop_size, new_overlap);
save_dir = sprintf('\\cropped%dX%dX%d', new_crop_size);
mkdir(filepath, save_dir);

crop_filepath = filepath;

if exist('preferences.mat', 'file')
    save('preferences.mat', 'crop_filepath', '-append');
else
    save('preferences.mat', 'crop_filepath');
end

step = new_crop_size - overlap_px;
abandoned_list = [];

for n = 1 : file_num
    file = strcat(filepath, file_name{n});
    img = imread3d(file);
    [height, width, depth] = size(img);
    
    if padding
        padding = new_crop_size;
        if height == new_crop_size(1)
            padding(1) = 0;
        end
        if width == new_crop_size(2)
            padding(2) = 0;
        end
        if depth == new_crop_size(3)
            padding(3) = 0;
        end
        img = pad_zero3d(img, padding, 'end');
    end
    
    
    
    % cut image
    cropped_num_current_stack = 0;
    
    idx_z = 0;
    for k = 1 : step(3) : (depth - floor(new_crop_size(3)/3) + 1)
        if k > depth-new_crop_size(3)+1
            if cover_all
                k = depth - new_crop_size(3) + 1;
            else
                break;
            end
        end
        idx_z = idx_z + 1;
        idx_y = 0;
        
        for i = 1 : step(1) : (height - floor(new_crop_size(1)/2) + 1)
            if i > height-new_crop_size(1)+1
                if cover_all
                    i = height - new_crop_size(1) + 1;
                    
                else
                    break;
                end
            end
            idx_y = idx_y + 1;
            
            idx_x = 0;
            for j = 1 : step(2) : (width - floor(new_crop_size(2)/2) + 1)
                if j > width-new_crop_size(2)+1
                    if cover_all
                        j = width - new_crop_size(2) + 1;
                    else
                        break;
                    end
                end
                idx_x = idx_x + 1;
                
                fprintf('processing img %d / %d : ', n, file_num)
                fprintf('block[%d  %d  %d]   ', idx_y, idx_x, idx_z);
                
                h = i + new_crop_size(1) - 1;
                w = j + new_crop_size(2) - 1;
                d = k + new_crop_size(3) - 1;
                fprintf('%d : %d, %d : %d, %d : %d\n', i, h ,j, w, k, d);
                block = img(i : h, j : w, k : d);
                
                block_id = sprintf('%03d-%06d.tif', n + block_bias, cropped_num_current_stack);
                pixel_sum = sum(block(:));
                pixel_var = var(double(block(:)));
                %block_sbr = sbr(block(:), 8, 0.2, 0.5);
                
                %writeSequence(block, [path save_dir sprintf('\\%d-%06d', n, cropped_num_current_stack)]);
                if ~save_all
                    if pixel_sum < pixel_threshold && pixel_var > var_threshold
                        fprintf('sum %d var %d : saved as %s', pixel_sum, pixel_var, block_id)
                        write3d(block, [filepath save_dir '\\' block_id], bitDepth);
                    else
                        abandoned_list = [abandoned_list; block_id];
                        fprintf('sum %d var %d : abandoned', pixel_sum, pixel_var)
                    end
                else  % save all blocks , no matter valid or not
                    fprintf('sum %d var %d : saved as %s', pixel_sum, pixel_var, block_id)
                    write3d(block, [filepath save_dir '\\' block_id], bitDepth);
                end
                fprintf('\n')
                cropped_num_current_stack = cropped_num_current_stack + 1;
            end
            
        end
        
    end
    cropped_num = cropped_num + cropped_num_current_stack;
end

%folder = [filepath save_dir];
%folder_new = [filepath save_dir sprintf('grid%d-%d-%d', idx_y, idx_x, idx_z)];
%dos(sprintf('ren %s %s', folder, folder_new));

if cropped_num == 0
    fprintf('valid blocks : 0 \nstack size smaller than crop size\n');
    return
end
grid_dim = [idx_y, idx_x, idx_z];
fprintf('valid blocks : %d / %d \n', cropped_num - size(abandoned_list, 1), cropped_num);

crop_size = new_crop_size;
overlap = new_overlap;
save('preferences.mat', 'crop_size', 'overlap', '-append');
if ~cover_all
    save('preferences.mat', 'grid_dim', '-append');
end
if ~save_all
    save('preferences.mat', 'abandoned_list','-append');
end

