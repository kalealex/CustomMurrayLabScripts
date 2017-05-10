function [ filteredArray, outlierIdx ] = AK_windowIQRfilt(dataArray, windowWidth, nIQRsOutlierCriterion )
%AK_windowIQRfilt is a custom outlier removal algorithm that removes and
%flags outliers in time series data, where an outlier is considered to be
%any value a given number * the IQR outside of the IQR, within a sliding
%temporal window.
%   INPUT:
%       dataArray: an array of time series data to be filtered
%       windowWidth: the number of elements in the array which represent
%           the width of a sliding temporal window; the data are filtered
%           based on the IQR of the subset of data within this window, so
%           the window should be broad enough not to filter out global
%           patterns in the data that represent signal
%       nIQRsOutlierCriterion: a criterion for the number of IQRs above or
%           below the 3rd and 1st quartiles, respectively, that a value at
%           the center of the sliding temporal window has to attain to be
%           considered an outlier in the dataset
%   OUTPUT:
%       filteredData: a filtered version of the dataArray where outliers
%           are replaced with nans
%       outlierIdx: a logical array indexing into dataArray, where 1=>
%           outlier and 0 => non-outlier

% setup output
filteredArray = dataArray;
outlierIdx = zeros(size(dataArray));

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
    % find 1st and 3rd quartile boundaries and IQR
    quartiles = quantile(dataArray(window),[.25 .75]);
    IQR = quartiles(2)-quartiles(1);
    % turn outliers into nans
    outlierIdx(window) = outlierIdx(window) + (dataArray(window) < (quartiles(1) - nIQRsOutlierCriterion*IQR) | dataArray(window) > (quartiles(2) + nIQRsOutlierCriterion*IQR));
    % advance position of temporal window
    iCenter = iCenter+1;
end

% apply filtering
outlierIdx = outlierIdx > 0;
filteredArray(outlierIdx) = nan;

end

