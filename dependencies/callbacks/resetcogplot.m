function resetcogplot(objh,~,scatterh)
% objh is the button handle
% scatterh is handle to associated scatter plot

if objh.Value
    len = length(scatterh.XData);
    scatterh.XData = NaN(len,1); % reset upon setting numofdots
    scatterh.YData = NaN(len,1); % reset upon setting numofdots
    scatterh.CData = linspace(0,1,len); % reset upon setting numofdots
end
objh.Value = 0; % reset

end

