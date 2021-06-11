function [image_out, pixels] = pad2OddSquare(image_in,plotting)
% this function resizes an input image to a square of odd-length (so we
% have a well defined center) by zero-padding

% initialize sz-vector and addline
sz = zeros(1,2);
addline = 0;

% get input image dimensions
[sz(1), sz(2)] = size(image_in); [~,maxindex] = max(sz);

% we need an extra line if the maxindex would result in an even square
if ~bitget(sz(maxindex),1) == 1 % true if maxindex = even number
    addline = 1;
end

% we pad the smaller dimension with zeroes to get a square image
if maxindex == 2
    image_out = padarray(image_in,[sz(2)-sz(1)+addline, addline],min(image_in(1,:)),'pre');
else
    image_out = padarray(image_in,[addline, sz(1)-sz(2)+addline],min(image_in(:,1)),'pre');
end

% output final length (len*len)
pixels = size(image_out,1);

% % optional plotting
if plotting == true
    figure;
    center = ceil(sz(maxindex)/2); % get center
    surf(image_out)
    shading interp
    view(2)
    axis image
    hold on
    plot3([center center],[1 sz(maxindex)],[1 1],'k','LineStyle','--','Linewidth',2)
    plot3([1 sz(maxindex)],[center center],[1 1],'k','LineStyle','--','Linewidth',2)
end

end

