function [output,lastpath] = GetFilePath(lastpath)
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

% filter_img = {'*.jpg; *.JPG; *.jpeg; *.JPEG; *.bmp; *.BMP; *.png; *.PNG; *.tif; *.TIF; *.tiff, *.TIFF',...
%     'Image Files (*.jpg,*.img,*.tiff,)'; ...
%     '*.jpg','jpg Files (*.jpg)';...
%     '*.JPG','JPG Files (*.JPG)';...
%     '*.jpeg','jpeg Files (*.jpeg)';...
%     '*.JPEG','JPEG Files (*.JPEG)';...
%     '*.bmp','bmp Files (*.bmp)';...
%     '*.BMP','BMP Files (*.BMP)';...
%     '*.png','png Files (*.png)';...
%     '*.PNG','PNG Files (*.PNG)';...
%     '*.tif','tif Files (*.tif)';...
%     '*.TIF','TIF Files (*.TIF)';...
%     '*.tiff','tiff Files (*.tiff)';...
%     '*.TIFF','TIFF Files (*.TIFF)'};
% 
% filter_video = {'*.avi; *.AVI; *.mp4; *.MP4; *.mpg; *.MPG; *.m4v; *.M4V; *.wmv; *.WMV',...
%     'Video Files (*.avi, *.mp4, *.mpg,)'; ...
%     '*.avi','avi Files (*.avi)';...
%     '*.AVI','AVI Files (*.AVI)';...
%     '*.mp4','mp4 Files (*.mp4)';...
%     '*.MP4','MP4 Files (*.MP4)';...
%     '*.mpg','mpg Files (*.mpg)';...
%     '*.MPG','MPG Files (*.MPG)';...
%     '*.m4v','m4v Files (*.m4v)';...
%     '*.M4V','M4V Files (*.M4V)';...
%     '*.wmv','wmv Files (*.wmv)';...
%     '*.WMV','WMV Files (*.WMV)'};
