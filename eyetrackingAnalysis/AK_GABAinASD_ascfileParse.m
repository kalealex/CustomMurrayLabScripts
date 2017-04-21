function [timepts,event] = AK_GABAinASD_ascfileParse(filename,startRow,endRow)
%AK_GABAinASD_ascfileParse parses asc files of eyelink data for further analysis
%   imports informtion from asc file: time points, X and Y pupil
%   coordinates, and pupil size
%   documents eyelink softeware judgements and descriptives of fixation for each time point
%   documents most recent event message at each time point

if ischar(filename) == 0; % check that filename is a string
    disp('filename is not a string')
end

% Import asc file
%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
%   column5: text (%s)
%	column6: text (%s)
%   column7: text (%s)
%	column8: text (%s)
%   column9: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%*s%[^\n\r]';

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

% check that data arrays are equal length, then cat
lDA = cellfun(@length,dataArray);
if range(lDA)==0 % if all same size
    % cat all
    data = [dataArray{1:end-1}];
else % if not same size
    for iDA = 1:length(dataArray)
        % pad each array
        clear pad
        pad = cell(max(lDA)-length(dataArray{iDA}),1); % pad rows with empty cells of proper size
        dataArray{iDA} = [dataArray{iDA};pad];
    end
    % then cat all
    data = [dataArray{1:end-1}];
end
    
%% remove extra information from data array

for iD = 1:length(data); % cycle through lines of data cell array
    if iscell(data(iD,1)) == 1; % check that input is a cell array
        if strcmp(data{iD,1},'MSG') && ~isempty(data{iD,1}) == 1;% search for MSG
             if regexp(data{iD,2},regexptranslate('wildcard','*RECORD*')) == 1 % search for RECORD
                 sRow = iD; % start row for data array at RECORD
             end
        elseif strcmp(data{iD,1},'END') && ~isempty(data{iD,1}) == 1
            eRow = iD; % end row for data array at END
        end
    end
end
% trim data array
try
    data = data(sRow:eRow,1:9); 
catch
    data = data(sRow:end,1:9);
end

%% Assign event and timepoints

% create output cell array: timepts
timepts = cell(length(data),6);
timepts(1,1:6) = {'time', 'X coordinate', 'Y coordinate', 'pupil size', 'fixation', 'event'};
tRow = 2; %set up rows for indexing into timepts

% create output structure of cell arrays: event
event.blink = cell(length(data),3); 
event.blink(1,1:3) = {'start time', 'end time', 'duration'};
event.fix = cell(length(data),6);
event.fix(1,1:6) = {'start time', 'end time', 'angular resolution', 'average X position', 'average Y position', 'average pupil size'};
event.sacc = cell(length(data),9);
event.sacc(1,1:9) = {'start time', 'end time', 'duration', 'start X position', 'start Y position', 'end X position', 'end Y position', 'amplitude (degrees)', 'velocity (degrees/second)'};
blinkRow = 2; % set up blink count variable 
fixRow = 2; % set up fixation count variable
saccRow = 2; % ste up accade count variable

for iL = 1:length(data); % cycle through lines of data cell array
    if isnumeric(str2double(data{iL,1})) && ~isnan(str2double(data{iL,1})) == 1 % search for numbers in the first column of data
        timepts(tRow,1) = (data(iL,1)); % store timepoint
        timepts(tRow,2) = (data(iL,2)); % store X coordinate
        timepts(tRow,3) = (data(iL,3)); % store Y coordinate
        timepts(tRow,4) = (data(iL,4)); % store pupil size
        tRow = tRow+1;
    elseif strcmp(data{iL,1},'MSG') == 1 % search for MSG
        clear S is_number
        S = data{iL,2}; % set S as string to examine
        is_number = zeros(1,length(S));
        for iS = 1:length(S); % cycle through positions within string
            is_number(iS) = isnumeric(str2double(S(iS))) && ~isnan(str2double(S(iS))) && str2double(S(iS))~=1i || isspace(S(iS)); % create index of (non 1i) numeric characters and spaces
        end
        timepts(tRow+1:end,6) = {S(is_number==0)}; % set event label equil to non-numaric string content
        timepts = timepts([1:tRow-1 tRow+1:end],1:6); % remove extra row, which marks event in data, from timepts array
    elseif regexp(data{iL,1},regexptranslate('wildcard','SBLINK R*')) == 1; % search for blinks
%         event.blink(blinkRow,1) = data(iL+1,1); % document time point of blink start
        timepts(tRow+1:end,5) = {nan}; % document blink
        timepts = timepts([1:tRow-1 tRow+1:end],:); % remove extra row, which marks event in data, from timepts array
    elseif regexp(data{iL,1},regexptranslate('wildcard','SFIX R*')) == 1 % search for fixations
%         event.fix(fixRow,1) = data(iL+1,1); % document time point of fixation start
        timepts(tRow+1:end,5) = {1}; % document fixation
        timepts = timepts([1:tRow-1 tRow+1:end],:); % remove extra row, which marks event in data, from timepts array
    elseif regexp(data{iL,1},regexptranslate('wildcard','SSACC R*')) == 1 % search for saccades
%         event.sacc(saccRow,1) = data(iL+1,1); % document time point of saccade start
        timepts(tRow+1:end,5) = {0}; % document saccade
        timepts = timepts([1:tRow-1 tRow+1:end],:); % remove extra row, which marks event in data, from timepts array
    elseif regexp(data{iL,1},regexptranslate('wildcard','EBLINK R*')) == 1 % search for blinks ending
        event.blink{blinkRow,1} = data{iL,1}(strfind(data{iL,1},'R ')+2:end); % document time point of blink start
        event.blink(blinkRow,2:3) = data(iL,2:3); % document blink end time point & duration
        blinkRow = blinkRow+1; % add one to blink count
        timepts = timepts([1:tRow-1 tRow+1:end],:); % remove extra row, which marks event in data, from timepts array
    elseif regexp(data{iL,1},regexptranslate('wildcard','EFIX R*')) == 1 % search for fixations ending
        event.fix{fixRow,1} = data{iL,1}(strfind(data{iL,1},'R ')+2:end); % document time point of fixation start
        event.fix(fixRow,2:6) = data(iL,2:6); % document fixation end time point, angular resolution, avg X, avg Y, & avg pupil
        fixRow = fixRow+1; % add one to blink count
        timepts = timepts([1:tRow-1 tRow+1:end],:); % remove extra row, which marks event in data, from timepts array
    elseif regexp(data{iL,1},regexptranslate('wildcard','ESACC R*')) == 1 % search for saccades ending
        event.sacc{saccRow,1} = data{iL,1}(strfind(data{iL,1},'R ')+2:end); % document time point of saccade start
        event.sacc(saccRow,2:9) = data(iL,2:9); % document saccade end time point, duration, start X, start Y, end X, end Y, amp, velocity
        saccRow = saccRow+1; % add one to blink count
        timepts = timepts([1:tRow-1 tRow+1:end],:); % remove extra row, which marks event in data, from timepts array
    end
end

% resize funtion output
tRow = tRow-1;
blinkRow = blinkRow-1;
fixRow = fixRow-1;
saccRow = saccRow-1;

timepts = timepts(1:tRow,1:6);
event.blink = event.blink(1:blinkRow,1:3);
event.fix = event.fix(1:fixRow,1:6);
event.sacc = event.sacc(1:saccRow,1:9);

end

