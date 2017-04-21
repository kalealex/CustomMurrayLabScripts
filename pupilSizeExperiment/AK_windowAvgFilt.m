function [ filteredArray ] = AK_windowAvgFilt(dataArray, windowWidth, avgFunction )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% setup output
filteredArray = dataArray;

iCenter = 1;
while iCenter < length(dataArray)
    % initialize temporal window
    if round(iCenter - windowWidth/2) < 1
        window = 1:round(iCenter + windowWidth/2);
    elseif round(iCenter + windowWidth/2) > length(dataArray)
        window = round(iCenter - windowWidth/2):length(dataArray);
    else
        window = round(iCenter - windowWidth/2):round(iCenter + windowWidth/2);
    end
    % center = avg inside of window
    filteredArray(iCenter) = avgFunction(dataArray(window));
    % advance position of temporal window
    iCenter = iCenter+1;
end

end

