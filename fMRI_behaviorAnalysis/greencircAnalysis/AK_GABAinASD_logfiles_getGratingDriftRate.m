function [dt, stimuli] = AK_GABAinASD_logfiles_getGratingDriftRate(filename)
%AK_GABAinASD_logfiles_getGratingDriftRate takes a log file name as input and returns
%the difference in time between successive stimulus events and an array a
%stimulus event names from the logfiile

if ischar(filename) == 0; % check that filename is a string
    disp('filename is not a string')
end

% import events and times from tab-delimited log file
%% Initialize variables.
delimiter = '\t';
startRow = 4;
endRow = inf;


%% Format string for each line of text:
%   column2: double (%f)
%	column3: text (%s)
%   column4: text (%s)
%	column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%s%s%f%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable

dataArray([1, 4]) = cellfun(@(x) num2cell(x), dataArray([1, 4]), 'UniformOutput', false);
data = [dataArray{1:end-1}];

% trim down to bmp names and time points where the event is a Picture
% object being displayed (rather than a Response or a Pulse)
data = data(strcmp(data(:,2),'Picture'),3:4);

% return vector of differences in time in milliseconds between successive
% stimulus events...
dt = (cell2mat(data(2:end,2)) - cell2mat(data(1:end-1,2))) ./ 10;
% and an array of stimulus event names
stimuli = data(:,1);

end

