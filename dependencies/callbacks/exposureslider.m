function exposureslider(hobj,~,CamObj)
% hobj is the button/slider handle
% CamObj is ueye cam object handle

if CamObj.CameraObj.IsOpened
    % get exposure range for these settings / camera
    exposure = CamObj.ExposureRange;
    
    % generate nonlinear vector
    range = [exposure.Minimum:exposure.Increment:exposure.Minimum*25,...
        linspace(exposure.Minimum*25,exposure.Maximum,75)];
    
    % set exposure
    CamObj.Exposure = range(round(hobj.Value));
else
    disp('Exposureslider cannot change exposure because cam is not open')
end
end

