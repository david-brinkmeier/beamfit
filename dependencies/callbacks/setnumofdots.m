function setnumofdots(objh,~,scatterh)
% objh is the button handle
% scatterh is handle to associated scatter plot

% convert textbox input to number and round
% if it is not a number 
if ~isnan(str2double(objh.String))
    len = round(str2double(objh.String));
    if len < 1
        len = 25;
        objh.String = '25';
        warning('Mininum number of dots is 1')
    end
    scatterh.XData = NaN(len,1); % reset upon setting numofdots
    scatterh.YData = NaN(len,1); % reset upon setting numofdots
    scatterh.CData = linspace(0,1,len); % reset upon setting numofdots
else
    warning('Input must be numerical, ignoring user input')
end

end
