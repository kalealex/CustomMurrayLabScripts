function [ filteredArray, outlierIdx ] = AK_windowIQRfilt(dataArray, windowWidth, nIQRsOutlierCriterion )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

