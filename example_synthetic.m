% Copyright (c) 2021 David Brinkmeier
% davidbrinkmeier@gmail.com
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.

clearvars
clc
close all

%% load dependencies
addpath(genpath('dependencies'))

%% Synthetic Example / Preparing input image

% Grid and pixel pitch
pixelpitch = 5.2; % in µm
pixels = 256;
range = [-pixelpitch*pixels/2, pixelpitch*pixels/2]; % in µm
xlimit = range(1):pixelpitch:range(2); % +/- µm
ylimit = range(1):pixelpitch:range(2); % +/- µm
[Xgrid,Ygrid] = meshgrid(xlimit,ylimit);

% use symmetric or elliptic gaussian beam
type = 'donut'; % beam model type can be 'symmetric' or 'elliptic' or 'donut'
fitvariant = 'donutgaussian'; % fit model can be 'gaussian' or 'donutgaussian' or 'donutgaussian-approx'

% beam specifications
w0 = 50; % beam radius (symmetric beam) in µm
w0x = 295; w0y = 125; % beam radius x/y (elliptic beam) in µm

% is the beam offset from the center (0,0) and rotated or donut?
x_offset = 0; % in µm
y_offset = -0; % in µm
rotation_angle = 40; % only applies to variant elliptic
Rtrepan = 150; % trepan radius in µm, only applies to donut gaussian 

%% Do we add noise? (most of the stuff below is only for filename when saving)
addnoise.enable = 1;
addnoise.saveimages = 0; % true enables saving images as png files
addnoise.quantum_well_depth = 2e1; % configures poisson / shot noise [recommended: 1e1-1e5] (high - low)
addnoise.sigma_read = 20; % configures readout noise / gaussian white noise [recommended: 0-50] (off - high)
addnoise.plot = true; % do we want to before / after noise?
addnoise.pixelpitch = pixelpitch;
if strcmp(type,'elliptic')
    addnoise.w0a = w0x;
    addnoise.w0b = w0y;
else
    addnoise.w0a = w0;
    addnoise.w0b = w0;
end
addnoise.x0 = x_offset;
addnoise.y0 = y_offset;
addnoise.length = pixels*pixelpitch; % size of image
addnoise.angle = rotation_angle;

% Calculate Beam Profile & add Noise if necessary
switch type
    case 'symmetric'
        % Variant 1: Symmetric gaussian, units in µm
        gaussian = 250*exp(-2*((sqrt((Xgrid-x_offset).^2+(Ygrid-y_offset).^2)).^2./w0^2));
        input_image = gaussian;
    case 'elliptic'
        % Variant 2: Elliptical gaussian, units in µm
        parameters = [250,x_offset,w0x,y_offset,w0y,-rotation_angle*pi/180];
        rotated_elliptical_gaussian = @(A,Xgrid,Ygrid) A(1)*exp( -2*(...
            ( Xgrid*cos(A(6))-Ygrid*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(A(3)^2) + ...
            ( Xgrid*sin(A(6))+Ygrid*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(A(5)^2) ) );
        gaussian_ellipse = rotated_elliptical_gaussian(parameters,Xgrid,Ygrid);
        input_image = gaussian_ellipse;
    case 'donut'
        % Trepanning Gaussian: Gaussian beam w0 on circular path (radius Rtrepan)
        % cf. Karkhin, Victor A., "Thermal Processes in Welding", 2019,
        % ISBN 978-981-13-5965-1, DOI 10.1007/978-981-13-5965-1,
        % cf. chapter 5.2.1.3 solution re-derived for Gaussian with
        % definition E ~ exp(-2(r/w)^2)
        donutgaussian = exp(-2*(((Xgrid-x_offset).^2+(Ygrid-y_offset).^2)+Rtrepan.^2)./w0^2).* ...
                            besseli(0,4.*sqrt((Xgrid-x_offset).^2+(Ygrid-y_offset).^2).*Rtrepan./w0^2,0);
        donutgaussian(isinf(donutgaussian) | isnan(donutgaussian)) = 0; % Replace NaNs and infinite values with zeros
        donutgaussian = 250.*donutgaussian./max(donutgaussian(:));
        if sum(donutgaussian(:)) == 0
            error('Bessel function of the first kind to avoid overflow or loss of accuracy, adjust ratio w0/Rtrepan to higher values')
        end
        input_image = donutgaussian;
    otherwise
        error('undefined type')
end

% add shot and readout noise
if addnoise.enable == 1
    input_image_noise = ShotAndReadNoise(input_image,addnoise,type);
    input_image_noise(input_image_noise > 250) = 250;
else
    input_image_noise = input_image;
end
input_img = input_image_noise;

%% Fit beam
fitmodels = {'gaussian','donutgaussian','donutgaussian-approx'};
if ismember(fitvariant,fitmodels)
    settings.fitvariant = fitvariant;
else
    error('Fitvariant has been set to unknown model %s, known models are: "%s", "%s", "%s".',fitvariant,fitmodels{:})
end

settings.CCDpixelpitch = pixelpitch; % Camera pixel pitch (square pixels) in µm
settings.scaleInputFactor = 1; % scale input image (value < 1 reduces pixel count -> faster computation, > 1 oversampling)
settings.useCOG = false; % true: use COG for all analysis, false: use index(max value) as start for regression, empty []: NONE
settings.Crop = 'off'; % Valid strings are 'off', 'user' and 'auto'
settings.AutoCropStrength = 0.1; % Accepted Values: 0-100 (scaled as 0% loss to 50% acceptable energy loss)
settings.PostProcessing = 'off'; % accepted values: 'off', 'low', 'medium', 'high', 'veryhigh', 'ultra', 'maximum'; Ultra adds kovesi, maximum adds TVDenoising
settings.ISO11146 = true; % Calculate 2nd order central moments and plot ellipse + output beam radius short/long
settings.plot = true; % Plot results / automatically true if we save figure
settings.view3D = true; % use 3D view instead of 2D (surf instead of imagesc), for video processing consider fixing zlim in guassfit.m
settings.zlim = []; % useful for plotting 3D view of videos, needed to avoid the surf plot from "jumping", ex. uint8 max = 255;
settings.colormap = 'jet'; % Matlab colormaps, i.e. 'gray', 'parula', 'jet'
settings.shading = 'flat'; % Shading: 'flat' or 'interp' recommended
settings.savefig = 0; % 0: off, 1: save PNG, 2: save PNG and FIG using optionally provided filename
settings.savedata = 0; % 0: off, 1: save MAT, 2: save MAT and CSV (CSV only contains only CSV-appropriate data)
settings.savebeamprofile = 0; % enables saving of all processed beam images inside result struct (disable to save memory)
settings.normalizedata = 1; % normalizes output maxval to maxval input
settings.visuals = 0; % 0: default, 1: no annotation, 2: no cross section, 3: no annotation or cross section
settings.filename = []; % If unused, initialize emtpy, i.e. as []
settings.outputpath = []; % if unused, initialize emtpy, i.e. as []
settings.progressbar = false; % enable or disable progressbars

results = beamfit(settings,input_img);

%% Check Toolbox dependencies
% core_functions = {'autocrop.m',...
%     'beamfit.m',...
%     'beamfitwrapper.m',...
%     'beamnoisefilt.m',...
%     'convert2grayscale.m',...
%     'draw_ellipse.m',...
%     'genORselectfigbyname.m',...
%     'GetFilePath.m',...
%     'image_moments.m',...
%     'noisecomp.m',...
%     'pad2OddSquare.m',...
%     'ShotAndReadNoise.m',...
%     'tight_subplot.m',...
%     'TVL1denoise.m',...
%     'UserCropImage.m'};
% 
% required_toolboxes = dependencies.toolboxDependencyAnalysis(core_functions);
% [fList,pList] = matlab.codetools.requiredFilesAndProducts('example_synthetic.m');
% license('inuse')

%% Deprecated

%     case 'donut'
%         % Approximation OK for ~ Rtrepan/w0 > 1.5
%         % Variant 3: This is a trepanning Gaussian beam as encountered in
%         % non-beam-rotating (non-beam-inverting) helical drilling /
%         % trepanning optics for helical drilling. A sufficiently fast
%         % trepanning beam w/ a long exposure camera looks results in the following img
%         donutgaussian_approx = 250*exp(-2*(sqrt((Xgrid-x_offset).^2+(Ygrid-y_offset).^2)-Rtrepan).^2./w0^2);
%         input_image = donutgaussian_approx;

% fprintf('w0 error: %3.2f %%\n',abs((1-w0/results.wx_radius))*100)
% fprintf('Rtrepan error: %3.2f %%\n',abs((1-Rtrepan/results.rtrepan))*100)