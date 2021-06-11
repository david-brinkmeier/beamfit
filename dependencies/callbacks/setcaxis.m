function setcaxis(objh,~,axh)
% objh is the button handle
% axh is axis handle which must be passed

% get handle to associated figure
% hfig = ancestor(objh, 'figure');  % or: gcbf

if objh.Value
    caxis(axh,[0 255])
else
    caxis(axh,'auto')
end

end

