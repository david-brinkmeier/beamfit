function results = beamfit(varargin)
% Copyright (c) 2021 David Brinkmeier
% davidbrinkmeier@gmail.com
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.

% [results,fighandles] = beamfit(varargin)
% This function is used to fit an elliptical gaussian beam or rotating symmetric gaussian to a picture
% taken by a CMOS/CCD Camera with square pixels

% Parse inputs
[settings, beam, background] = parse_inputs(varargin{:});
abort_flag = false;
maxval = 100; % normalize to 100 (%)

% Preallocate / clear output struct
results = struct();
% fighandles = struct(); & in case this code is ever rewritten properly

% Get screensize for plotting and positioning of figures / progressbars
screensize = get(0,'ScreenSize');

%% Prepare image (noise removal, resizing, zero padding, cropping)

% keep track of total calculation time
tic
if settings.BackgroundProvided
    beam = beam-background;
end

% If we donwsample then do it now for faster calculation
if settings.scaleInputFactor < 1 % resize image if scale value is not 1
    beam = imresize(beam,settings.scaleInputFactor); % resize image
end

% Crop image by user input if desired
% Autocrop will be used during Postprocessing at the appropriate time in beamnoisefilt.m
% UNLESS Postprocessing is disabled, then it will be used here
switch lower(settings.Crop) % convert string to lower case
    case 'off'
        % disp('cropping disabled')
    case 'user'
        beam = UserCropImage(beam,settings.colormap,true);
        if isempty(beam), return, end % when user aborts crop
    case 'auto'
        if strcmpi(settings.PostProcessing,'off')
            beam = autocrop(beam,settings.AutoCropStrength);
            warning('It is advised to use autocrop with at least PostProcessing = "low"!');
        end
    otherwise
        error('undefined crop type')
end

% Remove background noise
if settings.savebeamprofile, beam_noisy = beam; end % backup of original for output if we're saving it

if ~strcmpi(settings.PostProcessing,'off')
    [beam,abort_flag] = beamnoisefilt(beam,settings);
end

if abort_flag == true
    return
end

% If we use upsample then let's do it after noise removal and processing
if settings.scaleInputFactor > 1 % resize image if scale value is not 1
    beam = imresize(beam,settings.scaleInputFactor); % resize image
    if settings.savebeamprofile, beam_noisy = imresize(beam_noisy,settings.scaleInputFactor); end
end

% Pad to odd square
beam = pad2OddSquare(beam,false);
if settings.savebeamprofile, beam_noisy = pad2OddSquare(beam_noisy,false); end

% Calculate corrected pixel pitch dependent on Down/Upsampling
pixelpitch = settings.CCDpixelpitch/settings.scaleInputFactor; % interpolated pixel pitch

% normalize data
if settings.normalizedata
   beam = maxval.*beam./max(beam(:));
   if settings.savebeamprofile
      beam_noisy = maxval.*beam_noisy./max(beam_noisy(:));
   end
end

%% Prepare for fitting
% Prepare Grids for Fit
len = floor(max(size(beam))/2); len = len*pixelpitch;
[Xgrid,Ygrid] = meshgrid(-len:pixelpitch:len,-len:pixelpitch:len);

% calculate image moments
% im_moments = image_moments(beam,1); % only get center of gravity / first order moments
im_moments = image_moments(beam,0); % get center of gravity and central second order moments
% we need to interpolate the image_moments data to the correct values depending on pixel pitch and scale factor
im_moments.x = interp1(Xgrid(1,:),im_moments.x);
im_moments.y = interp1(Ygrid(:,1),im_moments.y);
im_moments.x1 = interp1(Xgrid(1,:),im_moments.x1);
im_moments.x2 = interp1(Xgrid(1,:),im_moments.x2);
im_moments.y1 = interp1(Ygrid(:,1),im_moments.y1);
im_moments.y2 = interp1(Ygrid(:,1),im_moments.y2);
im_moments.w = im_moments.w*pixelpitch; % minor
im_moments.l = im_moments.l*pixelpitch; % major 
im_moments.x_max = Xgrid(1,im_moments.x_max);
im_moments.y_max = Ygrid(im_moments.y_max,1);

if settings.useCOG == true % determine center of gravity and use this as starting value for regression
    % now we set input guess parameters for fitting
    Xstart = im_moments.x; % initial value x
    Ystart = im_moments.y; % initial value y
    Zstart = max(beam(:)); % initial value z / intensity
    DCstart = median(beam(:)); % initial value DC offset
elseif settings.useCOG == false % default to standard initial values
    Xstart = im_moments.x_max; % initial value x
    Ystart = im_moments.y_max; % initial value y
    Zstart = max(beam(:)); % initial value z / intensity
    DCstart = median(beam(:)); % initial value DC offset
else
    Xstart = 0; Ystart = 0; Zstart = 1; DCstart = 0;
end

%% ---Fit---
% Formula source: https://stackoverflow.com/questions/30420081/rotating-a-gaussian-function-matlab
% https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
% 2D Rotated Gaussian function (A requires 6 coefs). Input for 2d fit!

% ---Parameters---
bounds = round(max(Xgrid(:))/1); % we expect well-defined centered quadratic input (e.g. 3x3, 9x9, 301x301, etc.)
switch lower(settings.fitvariant)
    case 'gaussian'
        % fprintf('Xstart guess %3.1f µm, Ystart guess %3.1f µm\n',Xstart,Ystart)
        % Rotated gaussian fit function
        fitfun = @(A,X) A(1)*exp(-2*(...
                       (X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6))-A(2)*cos(A(6))+A(4)*sin(A(6))).^2/(A(3)^2) + ...
                       (X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6))-A(2)*sin(A(6))-A(4)*cos(A(6))).^2/(A(5)^2))) + ...
                        A(7);
        % inital (guess) parameters [Amp,x0,wx,y0,wy,theta,DC_offset]
        % wx/wy guess almost irrelevant because fit is very robust
        % A0 = [Zstart,Xstart,50,Ystart,50,0,DCstart];
        A0 = [Zstart,Xstart,50,Ystart,50,im_moments.theta,DCstart];
        % Define lower and upper bounds [Amp,xo,wx,yo,wy,fi]       
        % allowing for arbitrary angles ouside -pi/4 to pi/4 yields
        % consistently lower residuals but makes interpretation of X and Y
        % axis more ambigous
        lb = [0,-bounds,0,-bounds,0,-realmax('double'),0];
        ub = [realmax('double'),bounds,(bounds/2)^2,bounds,(bounds/2)^2,realmax('double'),max(beam(:))];
        % old variant
        % lb = [0,-bounds,0,-bounds,0,-pi/4,0];
        % ub = [realmax('double'),bounds,(bounds/2)^2,bounds,(bounds/2)^2,pi/4,max(beam(:))];
    
    case 'donutgaussian-approx'
        [Rtrepan,Xstart,Ystart] = get_Rguess(beam,pixelpitch,Xgrid(1,:),Ygrid(:,1));
        % fprintf('Rtrepan guess %3.1f µm, Xstart guess %3.1f µm, Ystart guess %3.1f µm\n',Rtrepan,Xstart,Ystart)
        % Trepanning rotationally symmetric gaussian fit function approx
        fitfun = @(A,X) A(1)*exp(-2*(sqrt((X(:,:,1)-A(2)).^2+(X(:,:,2)-A(3)).^2)-A(5)).^2./A(4)^2)+A(6);
        % inital (guess) parameters [Amp,x0,y0,w0,Rtrepan,DC_offset]
        A0 = [Zstart,Xstart,Ystart,50,Rtrepan,DCstart];
        % Define lower and upper bounds [Amp,x0,y0,w0,Rtrepan,DC_offset]
        lb = [0,-bounds,-bounds,0,0,0];
        ub = [realmax('double'),bounds,bounds,(bounds/2)^2,2*bounds,max(beam(:))];
        
    case 'donutgaussian'
        [Rtrepan,Xstart,Ystart] = get_Rguess(beam,pixelpitch,Xgrid(1,:),Ygrid(:,1));
        % fprintf('Rtrepan guess %3.1f µm, Xstart guess %3.1f µm, Ystart guess %3.1f µm\n',Rtrepan,Xstart,Ystart)
        % Trepanning rotationally symmetric gaussian fit function
        fitfun = @(A,X) A(6)+A(1).*exp(-2*(((X(:,:,1)-A(2)).^2+(X(:,:,2)-A(3)).^2)+A(5).^2)./A(4)^2).* ...
                                   besseli(0,4.*sqrt((X(:,:,1)-A(2)).^2+(X(:,:,2)-A(3)).^2).*A(5)./A(4)^2,0);
        % inital (guess) parameters [Amp,x0,y0,w0,Rtrepan,DC_offset]
        A0 = [Zstart,Xstart,Ystart,50,Rtrepan,DCstart];
        % Define lower and upper bounds [Amp,x0,y0,w0,Rtrepan,DC_offset]
        lb = [0,-bounds,-bounds,0,0,0];
        ub = [realmax('double'),bounds,bounds,(bounds/2)^2,2*bounds,max(beam(:))];
    
    otherwise
        error(['Fit model "',settings.fitvariant,'" unkown'])
end

% Numerical Grid for surface using lsqcurvefit requires wrapping X/Y meshgrid in 3D array
X(:,:,1) = Xgrid; X(:,:,2) = Ygrid;
options = optimset('Display','off'); % disable console output

if settings.progressbar == 1 % only show progressbar if requested
    progressbarfig = uifigure('Position',[round(0.5*(screensize(3)-400)) round(0.5*(screensize(4)-76)) 400 76]);
    progressbar = uiprogressdlg(progressbarfig,'Title','Computing Gaussian Fit','Indeterminate','on','Cancelable','on');
end

% Fit sample data -> execute lsqcurvefit
[A,resnorm,~] = lsqcurvefit(fitfun,A0,X,beam,lb,ub,options);
A(4)

% Close dialog box
if settings.progressbar == 1
    close(progressbar), close(progressbarfig), clear progressbar progressbarfig
end
tcalc = toc;

% debugging: evaluate fit and subtract beam
% figure; surf(fitfun(A,X)-beam), shading interp

%% Output data
fprintf(2,['Model:',32,'<strong>',upper(settings.fitvariant),'</strong>',32]);
results.filename = settings.filename;

if settings.savebeamprofile
    % save input image noisy / denoised and fitted function over same grid
    results.image = struct();
    results.image.filename = settings.filename;
    results.image.Xgrid = Xgrid;
    results.image.Ygrid = Ygrid;
    results.image.Z_denoised = beam;
    results.image.Z_input = beam_noisy;
    results.image.Fit = fitfun(A,X);
    if settings.BackgroundProvided, results.image.background = background; end
end

switch lower(settings.fitvariant)
    case 'gaussian'
        results.intensity = A(1);
        results.wx_pos = A(2);
        results.wy_pos = A(4);
        results.DC_offset = A(7);
        if abs(real(exp(1i*A(6)))) < abs(imag(exp(1i*A(6))))
            % if real part > imag part, then we must switch axes
            % this only really affects plotting of the projected lines
            % if the rotation real part exp^1i*angle is smaller than its
            % imaginary part what is defined as X-axis for the fit
            % is actually closer to the Y-axis of the camera reference
            % coordinate system -> then we switch the axes
            A(6) = A(6)+pi/2; % add 90 degree
            % now flip x and y
            results.wx_radius = A(5);
            results.wy_radius = A(3);
            A(3) = results.wx_radius;
            A(5) = results.wy_radius;
        else
            results.wx_radius = A(3);
            results.wy_radius = A(5);
        end
        A(6) = mod(A(6),pi);
        results.beam_angle = A(6);
        % this is just for plotting
        possible_angles = [180-A(6)*180/pi,90-A(6)*180/pi,-A(6)*180/pi];
        [~,idx_angle] = min(abs(possible_angles));
        fprintf('--- COG_x: %3.1f µm, COG_y: %3.1f µm, w0x: %3.1f µm, w0y %3.1f µm, tcalc = %2.2f s\n',...
            results.wx_pos,results.wy_pos,results.wx_radius,results.wy_radius,tcalc);
    case {'donutgaussian','donutgaussian-approx'}
        switch lower(settings.fitvariant)
            case 'donutgaussian'
                results.intensity = max(fitfun(A,X),[],'all');
                if A(5)/A(4) > 9
                  warning('Numerical solution of Bessel function results in numerical overflow and loss of accuracy, it is adviced to use donutgaussian-approx for Rtrepan/w0 > 9');  
                end
            case 'donutgaussian-approx'
                results.intensity = A(1);
                if A(5)/A(4) < 1.5
                    warning('approximate solution is inaccurate for Rtrepan/w0 < 1.5, use fit model "donutgaussian"!');
                end
        end
        results.wx_pos = A(2);
        results.wx_radius = A(4);
        results.wy_pos = A(3);
        results.wy_radius = A(4);
        results.rtrepan = A(5);
        results.DC_offset = A(6);
        fprintf('--- COG_x: %3.1f µm, COG_y: %3.1f µm, w0: %3.1f µm, Rtrepan %3.1f µm, tcalc = %2.2f s\n',...
            results.wx_pos,results.wy_pos,results.wx_radius,results.rtrepan,tcalc);
end

results.cogx = im_moments.x;
results.cogy = im_moments.y;
results.resnorm = resnorm; % resnorm = sum of squared residuals

if settings.ISO11146 == true
    results.w_short_11146 = im_moments.w;
    results.w_long_11146 = im_moments.l;
end

%% Plotting
% todo: Implementation acces using fighandles
% use fighandles to store handles to all surf/pcolor/line/annotations inside fitfig
% if figure exits, set data properly, i.e. https://stackoverflow.com/questions/13102654/how-should-i-update-the-data-of-a-plot-in-matlab
% right now for each subplot we clear the axes using cla
% and clear all annotations using delete(findall(fitfig,'type','annotation'))
% this is hacky, but at least it prevents memory leaks
% if we don't do this repeated plots to fitfig stack up and slow
% progressively slow down the execution
% making plot more efficient yields up to ~150% speed increase
% for 256px but only ~15% for 768px (for higher px count denoising and
% fitting of model has the biggest impact).
% thus, for now this is not a priority

if settings.plot || settings.savefig
    
    % ---Define Colors---
    if strcmpi(settings.colormap,'parula')
        color1 = 'k';
        color2 = '.r';
    elseif strcmpi(settings.colormap,'gray')
        color1 = 'g';
        color2 = '.r';
    elseif strcmpi(settings.colormap,'jet')
        color1 = 'm';
        color2 = '.w';
    end
    
    % ---Plot Data---
    [fitfig,figexists] = genORselectfigbyname('Beamfitting');
    % ex. access surface plot: fitfig.Children(3).Children
    % [fighandle.[3subplots],[surfaceandvariouslineplots]
    % directly write into this instead of redrawing..but lots of work
    
    if ~figexists
        set(fitfig,'Position',[screensize(3)-800-8 46 800 800]);
    end
    
    ax = subplot(4,4,[5,6,7,9,10,11,13,14,15]); cla
    if settings.view3D == true
        surf(Xgrid(1,:),Ygrid(:,1),beam);
        shading(gca, settings.shading); view(-14,50);
        zlabel('relative intensity')
    else
        pcolor(Xgrid(1,:),Ygrid(:,1),beam), shading(gca, settings.shading); %shading(gca, settings.shading);
        set(gca,'YDir','normal')
    end
    
    if settings.ISO11146 == true
        if settings.view3D == true % then we need to make grid and values available for draw_ellipse to interpolate z values
            im_moments.xgrid = Xgrid;
            im_moments.ygrid = Ygrid;
            im_moments.zvals = beam;
            im_moments.pixelpitch = pixelpitch;
        end
        hold on
        draw_ellipse(im_moments, 'LineStyle',':', 'color', color1, 'Linewidth', 2.0, 'elements', {'ellipse','major','minor'}) % Ellipsoid
        hold off
    end
    
    xlim([Xgrid(1) Xgrid(end)]), ylim([Ygrid(1) Ygrid(end)]);
    xlabel('distance in µm'), ylabel('distance in µm'), colormap(settings.colormap);
    box on, ax.LineWidth = 1; ax.FontSize = 12;
    
    % set ticks...uneven number of ticks to get the 0 tick
    tickvals = round(linspace(0,Xgrid(end)*0.95,4),-1);
    set(gca,'XTick',[fliplr(-tickvals(2:end)),tickvals])
    set(gca,'YTick',[fliplr(-tickvals(2:end)),tickvals])
    
    if ~isempty(settings.zlim)
        zlim([0 settings.zlim]);
    end
    
    if max(beam(:)) > 250
        title('warning, image may be overexposed!','Color','r','FontSize',14)
    end
    
    % set new aspect ratio
    AR = daspect;
    daspect([AR(1) AR(2) 2.5*AR(3)]);
    
    % Plot vertical and horizontal axis
    % obsolete (old implementation for -pi/4 to pi/4)
    % vx_h = Xgrid(1,:); vy_v=Ygrid(:,1);
    
    % generate points along _horizontal & _vertical axis
    interpvals = 255;
    switch lower(settings.fitvariant)
        case 'gaussian'
            % old implementation for -pi/4 to pi/4
            % problem: fit limited to -pi/4 to pi/4 because we can't handle tan(0)
            % vy_h = -tan(A(6))*(vx_h-A(2))+A(4); % gaussian y = ax+c
            % vx_v = -tan(A(6))*(A(4)-vy_v)+A(2); % gaussian y = ax+c

            % new implementation: generate vector with n samples between
            % the grid limits -len and len. rotate this vector which is
            % assumed to be horizontal [1,0] using exp(1i*angle*vector)
            % resulting real part is x values, imaginary part is y values
            % for the vertical line simply add +pi/2
            % this would be the horizontal line for abs(A(6)) < pi/4            
            horline = exp(1i*-A(6))*sqrt(2)*linspace(-2*len,2*len,interpvals)+(A(2)+1i*A(4));
            vx_h = real(horline);
            vy_h = imag(horline);
            % this would be the vertical line for abs(A(6)) < pi/4
            vertline = exp(1i*(-A(6)+pi/2))*sqrt(2)*linspace(-2*len,2*len,interpvals)+(A(2)+1i*A(4));
            vx_v = real(vertline);
            vy_v = imag(vertline);

        case {'donutgaussian','donutgaussian-approx'}
            % old implementation for -pi/4 to pi/4
            % vy_h = 0*vx_h+A(3); % donut gaussian
            % vx_v = 0*vy_v+A(2); % donut gaussian
            
            vx_h = linspace(-2*len,2*len,interpvals); 
            vy_v = linspace(-2*len,2*len,interpvals);
            vy_h = zeros(1,interpvals)+A(3); % donut gaussian
            vx_v = zeros(1,interpvals)+A(2); % donut gaussian
    end
    
    InterpMethod = 'linear'; % 'nearest','linear','spline','cubic'
    hPoints = interp2(Xgrid,Ygrid,beam,vx_h,vy_h,InterpMethod);
    vPoints = interp2(Xgrid,Ygrid,beam,vx_v,vy_v,InterpMethod);
    vcenter = interp2(Xgrid,Ygrid,beam,A(2),A(4));
    
    % plot lines
    if ismember(settings.visuals,[0, 1])
        hold on;
        if settings.view3D == true
            switch lower(settings.fitvariant)
                case 'gaussian'
                    plot3(A(2),A(4),vcenter,'+b'); % gaussian beam y-axis
                case {'donutgaussian','donutgaussian-approx'}
                    plot3(A(2),A(3),vcenter,'+b'); % donut gaussian x-axis
            end
            plot3(vx_h,vy_h,1.025.*hPoints,color2);
            plot3(vx_v,vy_v,1.025.*vPoints,color2);
        else
            switch lower(settings.fitvariant)
                case 'gaussian'
                    plot(A(2),A(4),'+b'); % gaussian beam
                case {'donutgaussian','donutgaussian-approx'}
                    plot(A(2),A(3),'+b'); % donut gaussian
            end
            plot(vx_h,vy_h,color2);
            plot(vx_v,vy_v,color2);
        end
        hold off;
    end
    
    % Plot cross sections
    dmin = 1.1*min(beam(:)); xfit = linspace(Xgrid(1,1),Xgrid(1,end),interpvals);
    dmax = 1.1*max(beam(:)); yfit = linspace(Ygrid(1,1),Ygrid(end,1),interpvals);
    
    switch lower(settings.fitvariant)
        case 'gaussian'
            hfit= A(7)+A(1)*exp(-2*(xfit-A(2)).^2/(A(3)^2)); % gaussian beam x-axis
            vfit= A(7)+A(1)*exp(-2*(yfit-A(4)).^2/(A(5)^2)); % gaussian beam y-axis
            xposh = (vx_h-A(2))/cos(A(6))+A(2); % rotated ellipse
            xposv = (vy_v-A(4))/cos(A(6))+A(4); % rotated ellipse
        case 'donutgaussian'                 
            hfit = A(6)+A(1).*exp(-2*(((xfit-A(2)).^2)+A(5).^2)./A(4)^2).* ...
                   besseli(0,4.*sqrt((xfit-A(2)).^2).*A(5)./A(4)^2,0); % donut gaussian x-axis
            vfit = A(6)+A(1).*exp(-2*(((yfit-A(3)).^2)+A(5).^2)./A(4)^2).* ...
                   besseli(0,4.*sqrt((yfit-A(3)).^2).*A(5)./A(4)^2,0); % donut gaussian y-axis
            xposh = vx_h; % donut
            xposv = vy_v; % donut
        case 'donutgaussian-approx'
            hfit = A(1)*exp(-2*(sqrt((xfit-A(2)).^2)-A(5)).^2./A(4)^2)+A(6); % donut gaussian x-axis approx
            vfit = A(1)*exp(-2*(sqrt((yfit-A(3)).^2)-A(5)).^2./A(4)^2)+A(6); % donut gaussian y-axis approx
            xposh = vx_h; % donut
            xposv = vy_v; % donut
    end

    % projection on horizontal axis
    axhor = subplot(4,4,[1,2,3]); cla;
    plot(xposh,hPoints,'r',xfit,hfit,'black','LineWidth',1.5); grid off; axis([-bounds,bounds,dmin,dmax]);
    axhor.Position = [0.1300 0.7 0.57 0.1566]; axis off
    
    % projection on vertical axis
    axver = subplot(4,4,[8,12,16]); cla;
    plot(vPoints,xposv,'r',vfit,yfit,'black','LineWidth',1.5); grid off; axis([dmin,dmax,-bounds,bounds]);
    axver.Position = [0.7075 0.122 0.1566 0.5700]; axis off
    
    if ismember(settings.visuals,[0, 2])
        % delete old annotations: not doing this causes a memory leak and poor performance
        delete(findall(fitfig,'type','annotation'))
        annot_offset = 0.2;
        % Annotation / Text Setup
        textsetup = {'FitBoxToText','off','EdgeColor','none','LineWidth',1,'FontName','Arial',...
            'BackgroundColor',[1 1 1],'Color',[0 0 0],'FaceAlpha',0,'FontSize',12,'FontWeight','normal'};
        
        switch lower(settings.fitvariant)
            case 'gaussian'
                str1 = ('I (a.u.) / x_0 [µm] / w_0x [µm] / y_0 [µm] / w_0y [µm] / angle [°]');
            case {'donutgaussian','donutgaussian-approx'}
                str1 = ('I (a.u.) / x_0 [µm] / w_0x [µm] / y_0 [µm] / w_0y [µm] / r_{trepan} [µm]');
        end
        annotation('textbox',...
            [0.13 0.645+annot_offset 0.56915 0.0625],...
            'String',str1,...
            'FontSize',12,...
            'FontName','Arial',...
            'FontWeight','normal',...
            'EdgeColor',[0 0 0],...
            'LineWidth',1,...
            'FaceAlpha',1.0,...
            'BackgroundColor',[1 1 1],...
            'Color',[0 0 0]);
        
        switch lower(settings.fitvariant)
            case 'gaussian'
                annotation('textbox',...
                    [0.131 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(1),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.21 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(2),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.29 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(3),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.385 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(4),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.47 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(5),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.57 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(possible_angles(idx_angle),'%.2f'),textsetup{:});
                    % 'String',num2str(-A(6)*180/pi,'%.2f'),textsetup{:});
                    % 'String',num2str(mod(-A(6)*180/pi,90),'%.2f'),textsetup{:});
            case {'donutgaussian','donutgaussian-approx'}
                annotation('textbox',...
                    [0.131 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str((results.intensity),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.21 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(2),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.29 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(4),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.385 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(3),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.47 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(4),'%.2f'),textsetup{:});
                annotation('textbox',...
                    [0.57 0.613+annot_offset 0.1 0.0625],...
                    'String',num2str(A(5),'%.2f'),textsetup{:});
        end
        
        if settings.ISO11146 == true
            str2 = strcat('ISO-11146: w =',32,'(',num2str(im_moments.l,'%.2f'),32,'/',...
                32,num2str(im_moments.w,'%.2f'),32,'µm), x_0',32,'/',32,'y_0 =',32,num2str(im_moments.x,'%.2f'),...
                32,'/',32,num2str(im_moments.y,'%.2f'),32,'µm');
            annotation('textbox',...
                [0.13 0.645+0.06225+annot_offset 0.56915 0.0315],...
                'String',str2,...
                'FitBoxToText','off',...
                'EdgeColor',[0 0 0],...
                'LineWidth',1,...
                'FontName','Arial',...
                'BackgroundColor',[1 1 1],...
                'Color',[0 0 0],...
                'FaceAlpha',1.0,...
                'FontSize',10,...
                'FontWeight','bold');
        end
    end
    drawnow
end

%% Save figure / data

if ~settings.forcenosave
    if any([settings.savefig settings.savedata])
        if ~isempty(settings.outputpath)
            output_folder = fullfile(settings.outputpath,'results_beamfit');
        else
            output_folder = fullfile(pwd,'results_beamfit');
        end
        if ~exist(output_folder,'dir')
            mkdir(output_folder)
            fprintf(['Saving results to location\n','%s','\n'],output_folder)
        end
        if ismember(settings.savefig,[1 2])
            if ~isempty(settings.filename)
                if settings.savefig == 1
                    saveas(fitfig,fullfile(output_folder,strcat(settings.filename,'.png')))
                end
                if settings.savefig == 2
                    saveas(fitfig,fullfile(output_folder,strcat(settings.filename,'.png')))
                    saveas(fitfig,fullfile(output_folder,strcat(settings.filename,'.fig')))
                end
            else
                if settings.savefig == 1
                    saveas(fitfig,fullfile(output_folder,'output.png'))
                end
                if settings.savefig == 2
                    saveas(fitfig,fullfile(output_folder,'output.png'))
                    saveas(fitfig,fullfile(output_folder,'output.fig'))
                end
            end
        end
        if settings.savedata
            if ~isempty(settings.filename)
                % output mat
                save(fullfile(output_folder,strcat(settings.filename,'.mat')),'results')
                if settings.savedata == 2
                    % then also output csv
                    if isfield(results,'image')
                        csv_output = rmfield(results,{'image';'filename'});
                    else
                        csv_output = results;
                    end
                    writetable(struct2table(csv_output), ...
                        fullfile(output_folder,strcat(settings.filename,'.csv')))
                end
            else
                save(fullfile(output_folder,'output.mat'),'results')
                if settings.savedata == 2
                    % then also output csv
                    if isfield(results,'image')
                        csv_output = rmfield(results,{'image';'filename'});
                    else
                        csv_output = rmfield(results,'filename');
                    end
                    writetable(struct2table(csv_output),fullfile(output_folder,'output.csv'))
                end
            end
        end
    end
end

end

%%  Function: parse_inputs
function [settings, beam, background] = parse_inputs(varargin)
narginchk(2,3)
settings = varargin{1};
% convert to grayscale
beam = convert2grayscale(varargin{2});
% remove NaN
if any(isnan(beam(:)))
    error('Image contains NaN.')
end

switch nargin
    case 2
        background = [];
        settings.BackgroundProvided = false;
    case 3
        background = convert2grayscale(varargin{3});
        settings.BackgroundProvided = true;
end

if ~isfield(settings, 'forcenosave')
    settings.forcenosave = false;
end

end

function [Rtrepan_guess,Xstart,Ystart] = get_Rguess(beam,pixelpitch,x_vect,y_vect)
% this function calulates a decend guess value for the trepan radius when
% fitting a donutgaussian or trepan gaussian beam

% extract all values above a threshold
extracted = double(beam > 0.5*max(beam(:))).*beam;
% this should be a noisy donut, calculate centers of gravity / first order moments
im_moments = image_moments(extracted,1);

% generate a grid for the square image and calculate radial distances from COG
[cols, rows] = meshgrid(1:length(beam), 1:length(beam));
radial_dist = (rows - round(im_moments.y)).^2 + (cols - round(im_moments.x)).^2;

% we want a radial logical ring mask with inner and outer radius, multiply
% this mask pointwise with the image and integrate. we expect the trepan radius of
% this donut to roughly coincide with the radial ringsection where this
% integral yields the maximum

% generate ring sections for evaluation
ring_sections = floor(linspace(1,length(beam),100));

energy = zeros(length(ring_sections)-1,1);
for i = 1:length(ring_sections)-1
    ringPixels  = radial_dist >= ring_sections(i).^2 & radial_dist <= ring_sections(i+1).^2;
    energy(i) = sum(extracted.*ringPixels,'all');
end
[~,max_idx] = max(energy);

% Output guess radius for trepan value as well as X and Y start values
Rtrepan_guess = pixelpitch*(ring_sections(max_idx+1)+ring_sections(max_idx))/2;
Xstart = x_vect(round(im_moments.x));
Ystart = y_vect(round(im_moments.y));

% code for debugging / visualization
% surf(double(ringPixels)), view(2), shading interp
% surf(double(extracted.*ringPixels)), view(3), shading interp
% surf(double(ringPixels)), view(2), shading interp
% pause(0.05)
end


%% develop old code

% figure;
% plot(vx_h,vy_h,'LineWidth',2); hold on
% plot(vx_v,vy_v,'LineWidth',2);
% % scatter(vx_h,vy_h); hold on
% % scatter(vx_v,vy_v);
% 
% % this would be the horizontal line for abs(A(6)) < pi/4
% testvect = exp(1i*(-A(6)))*sqrt(2)*linspace(-len,len,200)+(A(2)+1i*A(4));
% % logicalvect = (abs(real(testvect)) > len) | (abs(imag(testvect)) > len);
% % testvect = testvect(~logicalvect);
% 
% % this would be the vertical line for abs(A(6)) < pi/4
% testvect2 = exp(1i*(-A(6)+pi/2))*sqrt(2)*linspace(-len,len,200)+(A(2)+1i*A(4));
% % logicalvect2 = (abs(real(testvect2)) > len) | (abs(imag(testvect2)) > len);
% % testvect2 = testvect2(~logicalvect2);
% 
% % plot(real(testvect),imag(testvect))
% % plot(real(testvect2),imag(testvect2))
% scatter(real(testvect),imag(testvect))
% scatter(real(testvect2),imag(testvect2))
% 
% % testvectxy = [real(testvect); imag(testvect)];
% % testvectxy2 = [real(testvect2); imag(testvect2)];
