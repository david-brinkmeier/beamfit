function OUT = draw_ellipse(E, varargin)
%DRAW_ELLIPSE Draw ellipse
%   Ex. draw_ellipse(E, 'color', 'r', 'Linewidth', 1.5, 'elements', {'ellipse','major','minor'})
%   DRAW_ELLIPSE(E) draws the ellipse parametrized by the structure E.
%
%   IM.DRAW_ELLIPSE(..., 'n', N) draws the ellipse with N points. the
%   default value is N=100.
%
%   IM.DRAW_ELLIPSE(..., 'elements', ELM) draws specified elements. ELM can
%   be a string or a cell array of string containing the following values:
%       - 'ellipse' (default)
%       - 'major'
%       - 'minor'
%       - 'semimajor'
%       - 'direction'
%
%   DRAW_ELLIPSE(..., 'Property', P) applies pairs of properties to the
%   displayed ellipse. All properties of the plot group are applicable.
%
%   See also: IM.get_ellipse

% --- I/O handling ------------------------------------------------------

% --- Input parser

p = inputParser();
p.KeepUnmatched = true;
addRequired(p, 'E', @isstruct);
addOptional(p, 'n', 100, @isnumeric);
addOptional(p, 'elements', 'ellipse', @(x) ischar(x) || iscellstr(x));

parse(p, E, varargin{:});
in = p.Results;
un = [fieldnames(p.Unmatched) struct2cell(p.Unmatched)]';
un = un(:);

% --- Elements to display
if ischar(in.elements)
    in.elements = {in.elements};
end

% --- Prepare output
out = struct('ellipse', [], 'major', [], 'minor', [], 'direction', []);

% --- Computation & display -----------------------------------------------

% --- Ellipse

if ismember('ellipse', in.elements)
    
    % --- Prepare ellipse
    t = linspace(0, 2*pi, in.n);
    x = in.E.x + in.E.w*cos(t)*sin(in.E.theta)+in.E.l*sin(t)*cos(in.E.theta);
    y = in.E.y + in.E.l*sin(t)*sin(in.E.theta)-in.E.w*cos(t)*cos(in.E.theta);
    
    if isfield(E,'zvals') == 1
        z = interp2(E.xgrid,E.ygrid,E.zvals,x,y,'linear');
        out.ellipse = plot3(x, y, z.*1.025, un{:});
    else
        out.ellipse = plot(x, y, un{:});
    end
    
end

% --- Major axis

if ismember('major', in.elements)
    
    if isfield(E,'zvals') == 1
        x = linspace(-in.E.l*cos(in.E.theta)+in.E.x, in.E.l*cos(in.E.theta)+in.E.x, in.n);
        y = linspace(-in.E.l*sin(in.E.theta)+in.E.y, in.E.l*sin(in.E.theta)+in.E.y, in.n);
        z = interp2(E.xgrid,E.ygrid,E.zvals,x,y,'linear');
        plot3(x,y,z,un{:})
    else
        out.major = plot(in.E.x+[-1 1]*in.E.l*cos(in.E.theta), ...
            in.E.y+[-1 1]*in.E.l*sin(in.E.theta), un{:});
    end
    
elseif ismember('semimajor', in.elements)
    
    out.major = plot(in.E.x+[0 in.E.l*cos(in.E.theta)], ...
        in.E.y+[0 in.E.l*sin(in.E.theta)], un{:});
    
end

% --- Minor axis

if ismember('minor', in.elements)
    
    if isfield(E,'zvals') == 1
        x = linspace(in.E.x+in.E.w*sin(in.E.theta),in.E.x-in.E.w*sin(in.E.theta), in.n);
        y = linspace(in.E.y-in.E.w*cos(in.E.theta),in.E.y+in.E.w*cos(in.E.theta), in.n);
        z = interp2(E.xgrid,E.ygrid,E.zvals,x,y,'linear');
        plot3(x,y,z,un{:})
    else
        out.minor = plot(in.E.x+[-1 1]*in.E.w*sin(in.E.theta), ...
            in.E.y+[1 -1]*in.E.w*cos(in.E.theta), un{:});
    end
    
end

% --- Direction

if ismember('direction', in.elements)
    
    out.direction = scatter(in.E.x+in.E.l*cos(in.E.theta), ...
        in.E.y+in.E.l*sin(in.E.theta));
    
end


% --- Output --------------------------------------------------------------
if nargout
    OUT = out;
end