function [output_patch, output_patch_nofilt] = autocrop(input_image,target_loss,input_image_nofilt)
% Copyright (c) 2021 David Brinkmeier
% davidbrinkmeier@gmail.com
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.

% Function to automatically crop beam image based on center of gravity
% analysis and patch extraction considering energy loss

% Collect total image energy and center of gravity for cropping
total_energy = sum(input_image(:));
COG = image_moments(input_image,1); % get center of gravity / only calc 1st moments (barycenters)
COG = round([COG.y COG.x]);

% Maximum square size is determined by image boundaries
sz_limits = min([COG(1)-1, COG(2)-1, size(input_image,1)-COG(1)-1, size(input_image,2)-COG(2)-1]);

% Calculate results for the maximum viable patch first
max_output_patch = input_image(COG(1)-sz_limits:COG(1)+sz_limits,COG(2)-sz_limits:COG(2)+sz_limits);
max_patch_energy = sum(max_output_patch(:));
min_energy_loss = (1-max_patch_energy/total_energy)*100; % energy lost in %

if exist('input_image_nofilt','var')
    max_output_patch_nofilt = input_image_nofilt(COG(1)-sz_limits:COG(1)+sz_limits,COG(2)-sz_limits:COG(2)+sz_limits);
end

% target_loss defines / sets the energy loss threshold which is acceptable
% for auto croppping. We set maximum value to 50% of total energy
target_energy_loss = (0.5-(3/200)*min_energy_loss)*target_loss+min_energy_loss;

% If the energy loss is below our threshold, it is worth seeing if we can
% make the patch smaller!
% Initialize starting patch size and energy_loss for while loop
sz = 5;
energy_loss = 100; % initialize with 100%

while sz <= sz_limits && energy_loss > target_energy_loss
    output_patch = input_image(COG(1)-sz:COG(1)+sz,COG(2)-sz:COG(2)+sz);
    if exist('input_image_nofilt','var')
        output_patch_nofilt = input_image_nofilt(COG(1)-sz:COG(1)+sz,COG(2)-sz:COG(2)+sz);
    end
    patch_energy = sum(output_patch(:));
    energy_loss = (1-patch_energy/total_energy)*100; % energy lost in %
    sz = sz+4;
end

if exist('output_patch','var')
    disp(strcat('autocrop energy loss:',32,num2str(energy_loss-min_energy_loss),32,'%'));
else
    output_patch = max_output_patch;
    if exist('input_image_nofilt','var')
        output_patch_nofilt = max_output_patch_nofilt;
    end
    warning('Autocrop uses maximum width. Consider cropping manually.')
end

end

