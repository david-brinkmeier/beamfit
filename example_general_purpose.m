clc
close all

%% load dependencies
addpath(genpath('dependencies'))

%% Set up your settings for beamfitting here
% settings for beamfitting
settings.fitvariant = 'gaussian'; % Fit model: 'gaussian' or 'donutgaussian' (donutgaussian is a trepanning symmetric gaussian)
settings.CCDpixelpitch = 5.2; % Camera pixel pitch (square pixels) in µm
settings.scaleInputFactor = 1; % scale input image (value < 1 reduces pixel count -> faster computation, > 1 oversampling)
settings.useCOG = true; % true: use COG for all analysis, false: use index(max value) as start for regression, empty []: NONE
settings.Crop = 'user'; % Valid strings are 'off', 'user' and 'auto'
settings.AutoCropStrength = 0; % Accepted Values: 0-100 (scaled as 0% loss to 50% acceptable energy loss)
settings.PostProcessing = 'low'; % accepted values: 'off', 'low', 'medium', 'high', 'veryhigh', 'ultra', 'maximum'; Ultra adds kovesi, maximum adds TVDenoising
settings.ISO11146 = false; % Calculate 2nd order central moments and plot ellipse + output beam radius short/long
settings.plot = true; % Plot results / automatically true if we save figure
settings.view3D = true; % use 3D view instead of 2D (surf instead of imagesc), for video processing consider fixing zlim in guassfit.m
settings.zlim = []; % useful for plotting 3D view of videos, needed to avoid the surf plot from "jumping", ex. uint8 max = 255;
settings.colormap = 'parula'; % Matlab colormaps: 'gray', 'parula', 'jet'
settings.shading = 'flat'; % Shading: 'flat' or 'interp' recommended
settings.savefig = 0; % 0: off, 1: save PNG, 2: save PNG and FIG using optionally provided filename
settings.savedata = 0; % 0: off, 1: save MAT, 2: save MAT and CSV (CSV only contains only CSV-appropriate data)
settings.savebeamprofile = 1; % enables saving of all processed beam images inside result struct (disable to save memory)
settings.normalizedata = 1; % normalizes output maxval to 100
settings.visuals = 0; % 0: default, 1: no annotation, 2: no cross section, 3: no annotation or cross section
settings.filename = []; % If unused, initialize emtpy, i.e. as []
settings.outputpath = []; % if unused, initialize emtpy, i.e. as []
settings.progressbar = false; % enable or disable progressbars

% set up video processing (only used if you select a video or videos)
videosettings.startframe = 1; % first frame for processing
videosettings.sumframe = 1; % number of frames for time-averaging (1 means no averaging)
videosettings.lastframe = inf; % last video frame to process (set to inf for complete video)
videosettings.equalize = true; % if frames are averaged, equalize each frame value of each frame
videosettings.saveVideo = 0; % enable to save a video from processing
videosettings.FPS = 5; % video framerate

%% prompt user for files and execute beam fitting
if ~exist('lastpath','var'), lastpath = []; end
[files,lastpath] = GetFilePath(lastpath);
if ~isempty(files) % handle GetFilePath abort
    output = beamfitwrapper(files,settings,videosettings);
end