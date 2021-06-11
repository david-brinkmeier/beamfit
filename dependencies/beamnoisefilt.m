function [output_image,abort] = beamnoisefilt(input_image,settings)
% Copyright (c) 2021 David Brinkmeier
% davidbrinkmeier@gmail.com
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.

% Filter input image to prepare for gaussian beam fit
% if we don't have a background noise image we estimate noise dc offset by
% extracting 4 edge quadrants and looking for median value to be taken as
% background noise offset

% initialize TV Denoising setting
UseMedFilt = true;
abort = false;

% parse filter length
switch settings.PostProcessing
    case 'low'
        UseMedFilt = false;
    case 'medium'
        filtersize = 3;
    case 'high'
        filtersize = 5;
    case 'veryhigh'
        filtersize = 7;
    case 'ultra'
        % also adds kovesi
        filtersize = 5;
    case 'maximum'
        % also adds kovesi and TVD
        filtersize = 5;
end

%% Remove DC Offset
% Partition image into 16x16 Blocks, get median for each corner and
% assume the median of this to be the noise-DC offset of the image
[sz_y, sz_x] = size(input_image);
northwest = input_image(1:round(sz_y/4), 1:round(sz_x/4));
northeast = input_image(1:round(sz_y/4), end-round(sz_x/4)+1:end);
southwest = input_image(end-round(sz_y/4)+1:end, 1:round(sz_x/4));
southeast = input_image(end-round(sz_y/4)+1:end, end-round(sz_x/4)+1:end);
concatenated = [northwest northeast; southwest southeast];
dc_noise_offset = median(median(concatenated(:)));

if dc_noise_offset > 0.75*max(input_image(:))
    warning('Noise estimation suggests there is too much noise or not enough signal, disabling DC removal')
    dc_noise_offset = 0;
end

output_image = input_image-dc_noise_offset; % assume median value to be constant additive noise

%% Denoising steps

switch settings.PostProcessing
    case 'ultra'
        output_image = noisecomp(output_image, 2, 7, 2.5, 6, 1); output_image(output_image < 0) = 0;
        output_image = medfilt2(output_image,[filtersize filtersize]); % remove residual salt and pepper noise
        if strcmp(settings.Crop,'auto')
            output_image = autocrop(output_image,settings.AutoCropStrength);
        end
    case 'maximum'
        output_image = noisecomp(output_image, 2, 7, 2.5, 6, 1); output_image(output_image < 0) = 0;
        output_image = medfilt2(output_image,[filtersize filtersize]); % remove residual salt and pepper noise
        % TV Denoising if enabled
        if strcmp(settings.Crop,'auto')
            % we don't want to run TVDenoise on the uncropped image,
            % but it also shouldn't be pre-filtered. We run autocrop on
            % a filtered image but output an unfiltered patch of the
            % same size and position
            [~,output_image] = autocrop(medfilt2(output_image,[filtersize filtersize]),settings.AutoCropStrength,output_image);
        end
        [output_image, abort] = TVL1denoise(output_image, 0.625, 1000, settings.progressbar);
        output_image(output_image<0) = 0; % positivity constraint in case of outliers
        
    otherwise
        output_image(output_image<0) = 0; % positivity constraint in case of outliers
        if UseMedFilt == true
            output_image = medfilt2(output_image,[filtersize filtersize]); % remove residual salt and pepper noise
        end
        % DC offset is now removed, now we can auto-crop
        if strcmp(settings.Crop,'auto')
            output_image = autocrop(output_image,settings.AutoCropStrength);
        end
end

end