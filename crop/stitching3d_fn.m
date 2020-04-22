function stitching3d_fn(grid_dim, new_overlap, delete_source)
% % 3-D stitching
% params:
%   -grid_dim: 3-element list, number of blocks along [height, width, depth]
%   -overlap: 3-element list, percentage of overlaps along [height, width, depth]
%   -delete_source: boolean, whether to delete source blocks after stitching 

tic

path(path, '_utils')

%grid_dim   = [13, 11, 10];
%overlap    = [0.5, 0.5, 0.5];

if exist('preferences.mat', 'file')
    load preferences.mat
    if exist('stitch_filepath', 'var') && exist(stitch_filepath, 'dir')
        filepath = stitch_filepath;
    elseif exist('crop_filepath', 'var') 
        filepath = crop_filepath;
    else
        filepath = pwd;
    end
else
    filepath = pwd;
end

[files,  filepath] = uigetfile({ '*.tif'; '*.*'}, 'MultiSelect', 'on', 'Select blocks...', filepath);

if ~iscell(files)
    if files == 0
        return
    end
    files = {files};
end

stitch_filepath = filepath;
if exist('preferences.mat', 'file')
    save('preferences.mat', 'stitch_filepath', '-append');
else
    save('preferences.mat', 'stitch_filepath');
end


n_files = length(files);

sampleFile = fullfile(filepath, files{1});
sampleInfo = imfinfo(sampleFile);
bitdepth = sampleInfo.BitDepth
sample_block = imread3d(sampleFile);
block_size = size(sample_block);

if (length(block_size) == 2)
    block_size = [block_size 1];
end

overlap_px = floor( new_overlap .* block_size);
block_margin = overlap_px / 2;

stitching_size = block_size .* grid_dim - overlap_px .* (grid_dim - 1);
step = block_size - block_margin;

block_idx = 0;

n_volumes = ceil(n_files / prod(grid_dim));

for nf = 1 : n_volumes
    stitched = zeros(stitching_size);
    
    end_z = block_margin(3);
    for k = 1 : grid_dim(3) %% depth
        
        begin_z = end_z - block_margin(3) + 1;
        end_z =  begin_z + step(3) - 1;
        block_margin_z = block_margin(3);
        
        if (k == 1)
            end_z = block_size(3); %% margin blocks should not be trimmed
            block_margin_z = 0;
        end
        
        end_y = block_margin(1);
        for i = 1 : grid_dim(1) %% height
            
            begin_y = end_y - block_margin(1) + 1;
            end_y =  begin_y + step(1) - 1;
            block_margin_y = block_margin(1);
            
            if (i == 1)
                end_y = block_size(1);
                block_margin_y = 0;
            end
            
            end_x = block_margin(2);
            for j = 1 : grid_dim(2) %% width
                
                begin_x = end_x - block_margin(2) + 1;
                end_x = begin_x + step(2) - 1;
                block_margin_x = block_margin(2);
                if (j == 1 )
                    end_x = block_size(2);
                    block_margin_x = 0;
                end
                
                block_idx = block_idx + 1;
                if block_idx > n_files
                    break;
                end
                tmp_block = imread3d(fullfile(filepath, files{block_idx}));
                tmp_block = tmp_block(block_margin_y + 1 : end, block_margin_x + 1 : end, block_margin_z + 1 : end);
                stitched(begin_y : end_y, begin_x : end_x ,begin_z : end_z) = tmp_block;
                
                fprintf('fusing block %d  at region [%d : %d, %d : %d, %d : %d]\n', block_idx, begin_y, end_y, begin_x, end_x, begin_z, end_z);
                
            end
        end
    end
    
    stitched = stitched ./ max(stitched(:)) * (2^bitdepth - 1);
    
    %write3d(stitched, fullfile(filepath, sprintf('fused-%04d.tif', nf)), bitdepth);
    
    
    if (bitdepth == 16)
        stitched = uint16(stitched);
    else
        stitched = uint8(stitched);
    end
    slice_dir = [filepath sprintf('stack-%04d', nf)];
    mkdir(slice_dir);
    for d = 1 : size(stitched, 3)
        imwrite(stitched(:,:,d),  fullfile(slice_dir, sprintf('slice-%04d-.tif', d)));
    end
end
toc

if delete_source
    for i = 1 : length(files)
        delete(fullfile(filepath, files{i}));
    end
end
