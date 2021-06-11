function setcolormap(objh,~,axh)
% objh is the button handle
% axh is axis handle which must be passed
% get handle to associated figure
% hfig = ancestor(objh, 'figure');  % or: gcbf

% select string associated with value, convert to string and pass to
% colormap with handle axh
colormap(axh,char(objh.String(objh.Value)))
end

