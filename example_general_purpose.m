% Copyright (c) 2021 David Brinkmeier
% davidbrinkmeier@gmail.com
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.

clc
clearvars
close all

%% INFO
% Select one or multiple images OR select one or multiple videos for
% processing. Videos may be time-averaged.
%
% Video results yield a nested structure of per frame results 
% and results for the complete video.
% Especially if processing a complete video it is advised to disable
% settings.savebeamprofile in order not to bloat the results.mat
% Video processing only supports output of video and .mat file, if you need
% every fig as .fig and results as .csv provide single images

%% load dependencies
addpath(genpath('dependencies'))

%% Set up your settings for beamfitting here
% settings for beamfitting
settings.fitvariant = 'donutgaussian'; % Fit model: 'gaussian' or 'donutgaussian' (donutgaussian is a trepanning symmetric gaussian) or 'donutgaussian-approx' for Rtrepan/w0 >> 1.5
settings.CCDpixelpitch = 5.2; % Camera pixel pitch (square pixels) in µm
settings.scaleInputFactor = 1; % scale input image (value < 1 reduces pixel count -> faster computation, > 1 oversampling)
settings.useCOG = 1; % true: use COG as centroid for start of fit, false: use index(max value) as start for regression, empty []: NONE
settings.Crop = 'user'; % Valid strings are 'off', 'user' and 'auto'
settings.AutoCropStrength = 0; % Accepted Values: 0-100 (scaled as 0% loss to 50% acceptable energy loss)
settings.PostProcessing = 'low'; % accepted values: 'off', 'low', 'medium', 'high', 'veryhigh', 'ultra', 'maximum'; Ultra adds kovesi, maximum adds TVDenoising
settings.ISO11146 = 1; % Calculate 2nd order central moments and plot ellipse + output beam radius short/long
settings.plot = 1; % Plot results / automatically true if we save figure
settings.view3D = 0; % use 3D view instead of 2D (surf instead of imagesc), for video processing consider fixing zlim in guassfit.m
settings.zlim = []; % useful for plotting 3D view of videos, needed to avoid the surf plot from "jumping", ex. uint8 max = 255;
settings.colormap = 'gray'; % Matlab colormaps: 'gray', 'parula', 'jet'
settings.shading = 'flat'; % Shading: 'flat' or 'interp' recommended
settings.savefig = 2; % 0: off, 1: save PNG, 2: save PNG and FIG using optionally provided filename
settings.savedata = 0; % 0: off, 1: save MAT, 2: save MAT and CSV (CSV only contains only CSV-appropriate data)
settings.savebeamprofile = 1; % enables saving of all processed beam images inside result struct (disable to save memory)
settings.normalizedata = 1; % normalizes output maxval to 100
settings.visuals = 0; % 0: default, 1: no annotation, 2: no cross section, 3: no annotation or cross section
settings.filename = []; % If unused, initialize emtpy, i.e. as []
settings.outputpath = []; % if unused, initialize emtpy, i.e. as []
settings.progressbar = false; % enable or disable progressbars

% set up video processing (only used if you select a video or videos)
videosettings.startframe = 1; % first frame for processing
videosettings.sumframe = 500; % number of frames for time-averaging (1 means no averaging)
videosettings.lastframe = inf; % last video frame to process (set to inf for complete video)
videosettings.equalize = true; % if frames are averaged, equalize each frame value of each frame
videosettings.saveVideo = 1; % enable to save a video from processing
videosettings.compression = 1; % 0 is uncompressed AVI (LARGE!!!), 1 is mp4 @ max quality
videosettings.FPS = 5; % video framerate

%% prompt user for files and execute beam fitting
if ~exist('lastpath','var'), lastpath = []; end
[files,lastpath] = GetFilePath(lastpath);
if ~isempty(files) % handle GetFilePath abort
    output = beamfitwrapper(files,settings,videosettings);
end