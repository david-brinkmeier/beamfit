function [output,lastpath] = GetFilePath(lastpath)
% Copyright (c) 2021 David Brinkmeier
% davidbrinkmeier@gmail.com
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.

% provides optional last path so repeated calls open the last folder for selection
isvideo = false; % initialize

% filters for images and videos
filter_img = {'*.jpg; *.JPG; *.jpeg; *.JPEG; *.bmp; *.BMP; *.png; *.PNG; *.tif; *.TIF; *.tiff, *.TIFF',...
    'Image Files (*.jpg,*.img,*.tiff,)'}; ...
filter_video = {'*.avi; *.AVI; *.mp4; *.MP4; *.mpg; *.MPG; *.m4v; *.M4V; *.wmv; *.WMV',...
    'Video Files (*.avi, *.mp4, *.mpg,)'}; ...
filter = [filter_img; filter_video];

% prompt user to select image(s) or video(s)
[filename, pathname] = uigetfile(filter,'MultiSelect','on','Select a file / files',lastpath);

% handle case if user aborts
if isequal(filename,0)
    disp('No files selected / Load cancelled.')
    output = [];
    return
else
end

fullfilename = strcat(pathname,filename);
if iscell(filename)
    [~,filename,fileext] = cellfun(@fileparts,filename,'Un',false);
else
    [~,filename,fileext] = fileparts(filename);
end

% check if we have a video
if contains(filter_video{1},fileext)
    isvideo = true;
end

% transpose...easier to read / debug in workspace
if iscell(filename)
    output.fullfilename = fullfilename.';
    output.filename = filename.';
    output.fileext = fileext.';
else
    % output as cell
    output.fullfilename = {fullfilename};
    output.filename = {filename};
    output.fileext = {fileext};
end
output.path = pathname;
output.isvideo = isvideo;
lastpath = output.path;
end