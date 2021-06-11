% Copyright (c) 2021 David Brinkmeier
% davidbrinkmeier@gmail.com
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.

%% init
clc
% add ueyedotnet assembly on first run
if ~exist('ueyedotnet','var')
    ueyedotnet = NET.addAssembly('C:\Program Files\IDS\uEye\Develop\DotNet\signed\uEyeDotNet.dll');
end

%% load dependencies
addpath(genpath('dependencies'))

%% These settings may be changed but the important settings can be modified
% using the GUI so if you're not sure leave these as is

% settings for beamfitting which are not modified using UI
settings.useCOG = false; % true: use COG for all analysis, false: use index(max value) as start for regression, empty []: NONE
settings.Crop = 'off'; % Valid strings are 'off', 'user' and 'auto'
settings.AutoCropStrength = 0; % Accepted Values: 0-100 (scaled as 0% loss to 50% acceptable energy loss)
settings.ISO11146 = true; % Calculate 2nd order central moments and plot ellipse + output beam radius short/long
settings.zlim = 100; % useful for plotting 3D view of videos, needed to avoid the surf plot from "jumping", ex. uint8 max = 255;
settings.colormap = 'jet'; % Matlab colormaps, i.e. 'gray', 'parula', 'jet'
settings.shading = 'flat'; % Shading: 'flat' or 'interp' recommended
settings.savefig = 0; % 0: off, 1: save PNG, 2: save PNG and FIG using optionally provided filename
settings.savedata = 0; % 0: off, 1: save MAT, 2: save MAT and CSV (CSV only contains only CSV-appropriate data)
settings.savebeamprofile = 0; % enables saving of all processed beam images inside result struct (disable to save memory)
settings.normalizedata = 1; % normalizes output maxval to 100
settings.visuals = 0; % 0: default, 1: no annotation, 2: no cross section, 3: no annotation or cross section
settings.filename = []; % If unused, initialize emtpy, i.e. as []
settings.outputpath = []; % if unused, initialize emtpy, i.e. as []
settings.progressbar = false; % enable or disable progressbars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               DO NOT CHANGE ANYTHING BEYOND THIS POINT                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Connect to camera on first run
if ~exist('CamObj','var'), CamObj = uEyeObj; end

% Initialize camera
CamObj.Initialize;

% open imagesc and draw one frame
[ueyefig,ueyefigexists] = genORselectfigbyname('uEye Viewer');

% init limits for ROI-crop
limits = [];

if ~ueyefigexists
    axUEYE = axes(ueyefig); axUEYE.Position = [0.13 0.17 0.775 0.69];
    ueyefig.Position = [50,50,720,560]; ueyefig.Resize = 'off';
    im = imagesc(axUEYE,CamObj.GrabImage.'); axis off;
    ueyetitle = title(axUEYE,{'string1','string2'}); % generate multiline title
    % generate all required UI elements
    % toggle button to stop capturing / fitting
    StopExecute = uicontrol(ueyefig,'Style', 'togglebutton', ...
        'String', 'Stop', ...
        'Position', [10,10,80,30]);
    % toggle button for caxis
    uicontrol(ueyefig,'Style', 'togglebutton', ...
        'String', '0-255', ...
        'Position', [600,450,50,30], ...
        'Callback', {@setcaxis,axUEYE});
    % toggle button for colormap
    colormapselector = uicontrol(ueyefig,'style','popup','Position',[538,435,60,45],...
        'string',{'jet','parula','gray'},...
        'Callback', {@setcolormap,axUEYE});
    % slider for exposure control
    uicontrol('Parent',ueyefig,'Style','slider','Position',[95,70,558,20],...
        'value',1,'min',1,'max',100,'BackgroundColor',[1 1 1],...
        'Callback', {@exposureslider,CamObj});
    % make panel for selection of beamfitting settings
    beamfitpanel = uipanel('Parent',ueyefig,'Title','Beamfitting','FontSize',12,...
        'BackgroundColor','white','Position',[.16,.015,.7175,.1]);
    % toggle button for ROI choosing (disable after setting limits)
    chooseROI = uicontrol(beamfitpanel,'Style', 'togglebutton', ...
        'String','select ROI','Position', [5,6,80,24]);
    % dropdown for selection of fitvariant
    selectFitvariant = uicontrol(beamfitpanel,'style','popup','Position',[92,4,80,25],...
        'string',{'gaussian','donutgaussian'},'Value',1);
    % dropdown for selection of postprocessing strength
    selectPostprocessing = uicontrol(beamfitpanel,'style','popup','Position',[180,4,80,25],...
        'string',{'off','low','medium','high','veryhigh','ultra','maximum'},'Value',3);
    % dropdown for scaleFactor (convert to number)
    scaleFactor = uicontrol(beamfitpanel,'style','popup','Position',[268,4,80,25],...
        'string',{'1.0','0.75','0.5','0.33','0.25','0.15','0.1'},'Value',3);
    % toggle for selection of 3Dview vs 2D view (expected default)
    enable3DView = uicontrol(beamfitpanel,'Style', 'togglebutton', ...
        'String','3D','Position',[355,6,30,24]);
    % toggle for selection of 3Dview vs 2D view (expected default)
    enablePlot = uicontrol(beamfitpanel,'Style', 'togglebutton', ...
        'String','Plot','Position',[391,6,30,24],'Value',1);
    % toggle for execution of beamfitting
    startBeamfit = uicontrol(beamfitpanel,'Style', 'togglebutton', ...
        'String','F I T','Position',[427,6,80,24],'FontWeight','bold');
end

% get camera model
cam_model = char(CamObj.SensorInfo.SensorID);
cam_model = cam_model(cam_model ~= '_'); % remove underscores

% get/set pixel pitch
settings.CCDpixelpitch = double(CamObj.SensorInfo.PixelSize)/100; % in µm
while ~StopExecute.Value
    % Display frames until we press Stop
    try
        im.CData = CamObj.GrabImage.';
    catch
        disp('Image acquisition error (╯°□°）╯︵ ┻�?┻)')
    end
    ueyetitle.String = {sprintf('IDS %s, Pixelpitch %1.2f µm',cam_model,settings.CCDpixelpitch);...
        sprintf('Exposure [min-max]: %u-%u',min(im.CData(:)),max(im.CData(:)))};
    
    if chooseROI.Value
        [~,limits] = UserCropImage(im.CData,settings.colormap,true);
        chooseROI.Value = 0; % reset toggle button
    end
    
    if startBeamfit.Value
        % parse UI elements
        settings.view3D = enable3DView.Value;
        settings.plot = enablePlot.Value;
        settings.colormap = char(colormapselector.String(colormapselector.Value));
        settings.fitvariant = char(selectFitvariant.String(selectFitvariant.Value));
        settings.PostProcessing = char(selectPostprocessing.String(selectPostprocessing.Value));
        settings.scaleInputFactor = str2double(char(scaleFactor.String(scaleFactor.Value)));
        
        if ~isempty(limits)
            beam = UserCropImage(im.CData,[],[],limits);
        else
            beam = im.CData;
        end
        % execute beamfit and store results
        beamfitresults = beamfit(settings,beam);
        % update / make pointing viewer
        % short circuited...isgraphics not evaluated if first condition fails
        if ~exist('COGfig','var') || ~isgraphics(COGfig)
            COGfig = genORselectfigbyname('COG Viewer');
            COGfig.Position = [50,50,720,560];
            axCOG = axes(COGfig,'Position',[0.13 0.225 0.775 0.685]);
            scatterhandle = scatter(axCOG,NaN(25,1),NaN(25,1),100,linspace(0,1,25),'filled');
            ylabel(axCOG,'Position in µm'), xlabel(axCOG,'Position in µm'), COGtitle = title(axCOG,'init');
            axCOG.FontSize = 16; box(axCOG,'on'); axis(axCOG,'equal'); colormap(axCOG,flipud(gray));
            % generate all required UI elements
            COGpanel = uipanel('Parent',COGfig,'Title','Settings','FontSize',12,...
                'BackgroundColor','white','Units','pixels','Position',[92,10,306,50]);
            dotstext = uicontrol(COGpanel,'style','text','Position',[66 6 50 20],...
                'BackgroundColor','w','FontSize',10,'String','Dots:');
            numofdotsbox = uicontrol(COGpanel,'style','edit','Position',[110 5 50 24],...
                'String','25','Callback',{@setnumofdots,scatterhandle});
            useCOG = uicontrol(COGpanel,'Style','togglebutton', ...
                'String','useCOG','Position',[168 5 60 24]);
            uicontrol(COGpanel,'Style','togglebutton', ...
                'String','Reset','Position',[235 5 60 24],'Callback',{@resetcogplot,scatterhandle});
        end
        % update scatterplot circular buffer with new entries
        if ~useCOG.Value
            scatterhandle.XData = [scatterhandle.XData(2:end),beamfitresults.wx_pos];
            scatterhandle.YData = [scatterhandle.YData(2:end),beamfitresults.wy_pos];
        else
            scatterhandle.XData = [scatterhandle.XData(2:end),beamfitresults.cogx];
            scatterhandle.YData = [scatterhandle.YData(2:end),beamfitresults.cogy];
        end
        % set pointing viewer title
        COGtitle.String = getCOGtitle(scatterhandle.XData,scatterhandle.YData);
    end
    drawnow
end

try %#ok<TRYNC>
    StopExecute.Value = 0; % reset toggle button
    startBeamfit.Value = 0; % reset toggle button
    useCOG.Value = 0; % reset toggle button
end
% Close camera
CamObj.Close;