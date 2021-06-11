clear global
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
type = 'donut'; % type can be 'symmetric' or 'elliptic' or 'donut'
w0 = 64; % beam radius (symmetric beam) in µm
w0x = 195; w0y = 125; % beam radius x/y (elliptic beam) in µm

% is the beam offset from the center (0,0) and rotated or donut?
x_offset = 64; % in µm
y_offset = -102; % in µm
rotation_angle = 40; % only applies to variant elliptic
Rtrepan = 224; % trepan radius in µm, only applies to donut gaussian 

% Do we add noise?
addnoise.enable = 1;
addnoise.saveimages = 0; % true enables saving images as png files
addnoise.quantum_well_depth = 3e1; % configures poisson / shot noise [recommended: 1e1-1e5] (high - low)
addnoise.sigma_read = 30; % configures readout noise / gaussian white noise [recommended: 0-50] (off - high)
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

%% Calculate Beam Profile & add Noise if necessary
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
        % Variant 3: This is a trepanning Gaussian beam as encountered in
        % non-beam-rotating (non-beam-inverting) helical drilling /
        % trepanning optics for helical drilling. A sufficiently fast
        % trepanning beam w/ a long exposure camera looks results in the following img
        donutgaussian = 250*exp(-2*(sqrt((Xgrid-x_offset).^2+(Ygrid-y_offset).^2)-Rtrepan).^2./w0^2);
        input_image = donutgaussian;
    otherwise
        error('undefined type')
end

% add shot and readout noise
if addnoise.enable == 1
    input_image_noise = ShotAndReadNoise(input_image,addnoise);
    input_image_noise(input_image_noise > 250) = 250;
else
    input_image_noise = input_image;
end
input_img = input_image_noise;

%% Fit beam

if strcmpi(type,'donut')
    settings.fitvariant = 'donutgaussian'; % Fit model: 'gaussian' or 'donutgaussian'
else
    settings.fitvariant = 'gaussian'; % Fit model: 'gaussian' or 'donutgaussian'
end

settings.CCDpixelpitch = pixelpitch; % Camera pixel pitch (square pixels) in µm
settings.scaleInputFactor = 1; % scale input image (value < 1 reduces pixel count -> faster computation, > 1 oversampling)
settings.useCOG = false; % true: use COG for all analysis, false: use index(max value) as start for regression, empty []: NONE
settings.Crop = 'user'; % Valid strings are 'off', 'user' and 'auto'
settings.AutoCropStrength = 0; % Accepted Values: 0-100 (scaled as 0% loss to 50% acceptable energy loss)
settings.PostProcessing = 'high'; % accepted values: 'off', 'low', 'medium', 'high', 'veryhigh', 'ultra', 'maximum'; Ultra adds kovesi, maximum adds TVDenoising
settings.ISO11146 = true; % Calculate 2nd order central moments and plot ellipse + output beam radius short/long
settings.plot = true; % Plot results / automatically true if we save figure
settings.view3D = false; % use 3D view instead of 2D (surf instead of imagesc), for video processing consider fixing zlim in guassfit.m
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

for i = 1:1
% Process image and save results
% input_img = input_image_noise+50*rand(size(input_image_noise));
input_img = input_image_noise;
% input_img = UserCropImage(extra_noisy_img,'gray',true,limits);
results = beamfit(settings,input_img); % gaussfit(settings,input_image_noise,background);
end

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