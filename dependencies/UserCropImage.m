function [image_out, limits] = UserCropImage(varargin)
% Prompts user to draw a rectangle, region enclosed by rectangle will be
% extracted from image

% input arguments [image,colormap,forcesquare,limits]
% if limits are provided, then we don't prompt user and simply use these
% limits to extract the image
% limits must be a struct where limits.x and limits.y are integer (2,1) arrays

% ex. get limits from user
% [~,limits] = UserCropImage(inputimage,colormap,true);

% Parse inputs
[input_image, cmap, ForceSquare, limits, prompt] = parse_inputs(varargin{:});
[max_y,max_x] = size(input_image);

if prompt
    cropfig = figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(input_image), axis image; set(gca,'Ydir','normal'), colormap(cmap), axis off    
    title('Instructions: Draw square and double click in square when finished!',...
        'FontSize',16,'Color','r','FontWeight','bold')
    
    if ForceSquare
        h = imrect('PositionConstraintFcn', @(x) [x(1) x(2) min(x(3),x(4))*[1 1]]); %#ok<IMRECT>
    else
        h = imrect;  %#ok<IMRECT>
    end
    
    rect = wait(h);
    if ~isempty(rect)
        limits.x = [round(rect(1)), round(rect(1)+rect(3))];
        limits.y = [round(rect(2)), round(rect(2)+rect(4))];
        if limits.x(1) < 1, limits.x(1) = 1; end
        if limits.x(2) > max_x, limits.x(2) = max_x; end
        if limits.y(1) < 1, limits.y(1) = 1; end
        if limits.y(2) > max_y, limits.y(2) = max_y; end
        close(cropfig)
    else
        warning('User aborted cropping')
        image_out = []; limits = [];
        return
    end
end




if ~isempty(limits)
    image_out = input_image(limits.y(1):limits.y(2),limits.x(1):limits.x(2),:);
end

end

function [input_image, cmap, ForceSquare, limits, prompt] = parse_inputs(varargin)
narginchk(3,4)

% input arguments [image,colormap,forcesquare,limits]
% if limits are provided, then we don't prompt user and simply use these
% limits to extract the image
% limits must be a struct where limits.x and limits.y are integer (2,1) arrays

% input image is argin 1
input_image = varargin{1};

% colormap is argin 2
cmap = varargin{2};

switch nargin
    case 3
        prompt = true;
        ForceSquare = varargin{3};
        limits.x = 0;
        limits.y = 0;
    case 4
        prompt = false;
        ForceSquare = varargin{3};
        limits = varargin{4};
end
end