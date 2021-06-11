function output = beamfitwrapper(files,settings,videosettings)
% this is a helper function which automates usage of beamfit
% it expects as input the regular settings struct for beamfit
% and an additional struct for video settings which is required if videos
% are processed
%
% set up video processing (only used if you select a video or videos)
% videosettings.startframe = 1; % first frame for processing
% videosettings.sumframe = 60; % number of frames for time-averaging (1 means no averaging)
% videosettings.lastframe = 120; % last video frame to process (set to inf for complete video)
% videosettings.equalize = true; % if frames are averaged, equalize each frame value of each frame
% videosettings.saveVideo = 0; % enable to save a video from processing
% videosettings.FPS = 5; % video framerate
%
%%

if ~files.isvideo
    % initialize structure to store results for all images
    if settings.normalizedata, settings.zlim = 100; end
    output = struct(); output.settings = settings;    
    for i = 1:length(files.fullfilename)
        % input is images
        fprintf('Processing image %g of %g\n',i,length(files.fullfilename))
        settings.outputpath = files.path; % if unised, initialize emtpy, i.e. as []
        settings.filename = files.filename{i}; % If unused, initialize emtpy, i.e. as []
        output.(sprintf('im_%g',i)) = beamfit(settings,imread(files.fullfilename{i}));
    end
else
    % input is videos
    settings.zlim = 100; % we normalize to 100 in getnextframe so we set zlimit to 100
    % initialize structure to store results for all images
    output = struct(); output.settings = settings;
    % if user cropping is requested we disable for beamfit and ask user once per video
    if strcmpi(settings.Crop,'user')
        settings.Crop = 'off';
        askCropOnce = true;
        skipCrop = false;
    else
        askCropOnce = false;
        skipCrop = false;
    end
    for i = 1:length(files.fullfilename)
        fprintf('Processing video %g of %g\n',i,length(files.fullfilename))
        % set filename / output folder
        settings.filename = [files.filename{i},files.fileext{i}];
        settings.outputpath = files.path;
        % disable saving because we don't export results for each frame
        % instead we save for the whole video afterwards
        settings.forcenosave = 1;
        % initialize structure to store results for current video
        vidresults = struct(); 
        vidresults.settings = settings; vidresults.videosettings = videosettings;
        % init video read
        vid = VideoReader(files.fullfilename{i}); %#ok<TNMLP>
        vidresults.videoproperties = get(vid); % save vid properties
        % get number of frames..fairly recent addition to VideoReader so we
        % check if we need to calculate it ourselves
        if isprop(vid,'NumFrames')
            framecount = vid.NumFrames;
        else
            framecount = floor(vid.Duration*vid.FrameRate);
        end
        % set lastframe to framecount if it exceeds framecount
        if videosettings.lastframe > framecount, videosettings.lastframe = framecount; end
        if videosettings.lastframe == videosettings.startframe, videosettings.lastframe = videosettings.startframe+1; end
        % frameindexes defines the indexes between which we temporally integrate
        frameindexes = (videosettings.startframe:videosettings.sumframe:videosettings.lastframe);
        % get a frame at a random position in the video and have user set crop limits based on this image
        % we use random position in case user starts repeatedly to get an idea of the beams range of motion
        if askCropOnce && ~skipCrop
            [~,limits] = UserCropImage(getnextframe(vid,frameindexes,videosettings.equalize,...
                floor(1+(length(frameindexes)-1)*rand(1))),settings.colormap,true);
            if isempty(limits), return, end % handle user abort crop
            if i == 1 && length(files.fullfilename) > 1
                % if the first video is processed and there are multiple videos to process
                % ask user if we want to use these crop limits for all subsequent videos
                skipCrop = askCropLimits();
            end
        end
        % prepare video if user wants to export video
        if videosettings.saveVideo
            vidwriterObj = VideoWriter([files.path, files.filename{i},'_beamfit','.mp4'],'MPEG-4'); %#ok<TNMLP>
            vidwriterObj.Quality = 100; % Quality setting.
            vidwriterObj.FrameRate = videosettings.FPS; % How many frames per second.
            open(vidwriterObj);
        end
        % integrate image over frame sections, crop and fit
        % we could crop first but not without some rewriting and
        % computational cost is practically irrelevant
        for j = 1:length(frameindexes)-1
            fprintf('Video %g/%g: Frame %g/%g\n',i,length(files.fullfilename),j,length(frameindexes)-1)
            % get next frame and crop if limits exist
            [beam,indicator] = getnextframe(vid,frameindexes,videosettings.equalize,j);
            % crop based on known limits if provided
            if askCropOnce, beam = UserCropImage(beam,[],[],limits); end
            % fit beam
            vidresults.frameresults.(indicator) = beamfit(settings,beam);
            % write current frame if we export video
            if videosettings.saveVideo, writeVideo(vidwriterObj, getframe(gcf)); end
            if j == 1 % if it is the first analysis of that video, check which fields exist
                % get all field names in struct
                allfields = fieldnames(vidresults.frameresults.(indicator));
                % remove fields which are 'filename' or 'image'
                allfields = allfields(~ismember(allfields,{'filename','image'}));
            end
            % write all scalars as vectors to output
            for k = 1:length(allfields)
                vidresults.(allfields{k})(j) = vidresults.frameresults.(indicator).(allfields{k});
            end
        end
        % close video if we export video
        if videosettings.saveVideo, close(vidwriterObj); end
        % save / export results from current video
        output.(sprintf('vid_%g',i)) = vidresults;
        if ismember(settings.savedata,[1 2])
            save([files.path, files.filename{i},'_beamfit','.mat'],'vidresults')
        end
    end
    % reset forcenosave
    settings.forcenosave = 0;
end
fprintf('Done!\n')

end

function skipCrop = askCropLimits()
answer = questdlg('Use these crop limits for all other videos?', ...
    'Crop settings', ...
    'Yes','No','No'); % choice,choice,defaultchoice
switch answer
    case 'Yes'
        skipCrop = true;
    case 'No'
        skipCrop = false;
end
end

function [output,framestring] = getnextframe(vid,frameindexes,equalize,index)
% this function takes an object vid (type VideoReader) and outputs
% time-averaged frames [frameindexes(1) through frameindexes(2)].
% frameindexes is defined as frameindexes = (startframe:sumframe:lastframe);
% index selects the appropriate range from frameindexes
% maximum value for index is length(frameindexes);
% output data type id double

getframeidx = @(vect,idx) [vect(idx) vect(idx+1)-1];
% this function integrates a 4D stack of RGB frames along the temporal dimension
% normalizes to 100 and outputs the resulting image as grayscale double
% single frames pass directly and get converted to grayscale double
current_indexes = getframeidx(frameindexes,index);
frames = read(vid,current_indexes);
if max(frames(:)) > 250
%     warning('frames in this video are overexposed')
end
if size(frames,4) > 1
    output = zeros(size(frames,1:2));
    for i = 1:size(frames,4)
        if equalize
            currframe = convert2grayscale(frames(:,:,:,i));
            currframe = 100*currframe./max(currframe(:));
            output = output+currframe;
        else
            output = output+convert2grayscale(frames(:,:,:,i));
        end
    end
    output = 100*output./max(output(:));
else
    output = convert2grayscale(frames);
    output = 100*output./max(output(:));
end
framestring = sprintf('frame_%u_%u',current_indexes(1),current_indexes(2));

% fprintf('resulting image is time-averaged between frame %3.0f and %3.0f\n',current_indexes(1),current_indexes(2))
% imagesc(output)
end