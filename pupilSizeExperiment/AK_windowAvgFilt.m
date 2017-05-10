function [ filteredArray ] = AK_windowAvgFilt(dataArray, windowWidth, avgFunction )
%AK_windowAvgFilt is a smoothing filter for time series data that replaces
%each data point at the center of a sliding temporal window with the
%average of the points within the window.
%   INPUT:
%       dataArray: an array of time series data to be filtered
%       windowWidth: the number of elements in the array which represent
%           the width of a sliding temporal window; the data are filtered
%           based on the average of the subset of data within this window,
%           so the larger the window is, the more dramatic the smoothing
%           from the filter will be; the window should be narrow enough to
%           leave the signal intact but broad enough not to be trivial
%       avgFunction: this is the function that will be used for averaging
%           within the temporal window; either @nanmean or @nanmedian are
%           recommended so that a single nan in dataArray doesn't turn the
%           filtered data into many nans
%   OUTPUT:
%       filteredArray: a copy of dataArray which has been smoothed using
%           the window average filter 

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

