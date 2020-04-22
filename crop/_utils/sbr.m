%% Calculate the signal-to-background ratio of an image.
%  SBR is defined as the ratio of signals(the mean of the top 20% pixels) to the
%  background(the mean of the bottom 50% pixels)
%  
%% Params:
%  -perc_sgl : percentage of signal pixels in the total pixel number.
%  -perc_bg  : percentage of background pixels in the total pixel number.

%% Expamles:
% image = imread3d('test.tif');
% bitdepth = 8;
% perc_sgl = 0.2;
% perc_bg = 0.7; % for light sheet fluorenscent mouse brain images 

%%
function I = sbr(image, bitdepth, perc_sgl ,perc_bg)
    assert (bitdepth == 8 || bitdepth == 16);
    bins = 2 ^ bitdepth;
    img_histo = zeros(1, bins);
    
    [rols, cols, depth] = size(image);
    n_pixels = rols * cols * depth;
    
    image = uint8(image);
    for d = 1 : depth
        for r = 1 : rols
            for c = 1 : cols
                img_histo(image(r,c,d) + 1) = img_histo(image(r,c,d) + 1) + 1;
            end
        end
    end
    
    % calcucate the singals
    cursor = bins;
    tmp = 0; % tmp is the number of top 0.5% brightest pixels, which is excluded from the signals 
    while tmp < n_pixels * 0.05
        tmp = tmp + img_histo(cursor);
        cursor = cursor - 1;       
    end
    
    n_signals = 0;% nums of signal pixels 
    signal_intensity = 0;
    while n_signals + tmp < n_pixels * (perc_sgl + 0.05)
        n_signals = n_signals + img_histo(cursor);
        signal_intensity = signal_intensity + cursor * img_histo(cursor);
        cursor = cursor - 1;       
    end
    signal_intensity = signal_intensity / n_signals;
    
    % calculate the background 
    background_intensity = 0;
    cursor = 1;
    n_bgs = 0;
    while n_bgs < n_pixels * perc_bg
        n_bgs = n_bgs + img_histo(cursor);
        background_intensity = background_intensity + cursor * img_histo(cursor);
        cursor = cursor + 1;
    end
    background_intensity = background_intensity / n_bgs;
    
    I = signal_intensity / background_intensity;
    
        
            
        

