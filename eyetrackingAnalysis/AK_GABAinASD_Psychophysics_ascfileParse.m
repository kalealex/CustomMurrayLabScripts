function [block] = AK_GABAinASD_Psychophysics_ascfileParse(filename,startRow,endRow)
%AK_GABAinASD_ascfileParse parses asc files of eyelink data for further analysis
%   imports informtion from asc file: time points, X and Y pupil
%   coordinates, and pupil size
%   partitions data array into segments representing blocks of experiment
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

%% Create original data array

Odata = [dataArray{1:end-1}]; % original data

%% remove extra information from data array and create output structure

sCount = 0; % create start count variable
eCount = 0; % create end count variable

for iD = 1:length(Odata); % cycle through lines of data cell array
    if iscell(Odata(iD,1)) == 1; % check that input is a cell array
        if strcmp(Odata{iD,1},'MSG') && ~isempty(Odata{iD,1}) == 1;% search for MSG
             if regexp(Odata{iD,2},regexptranslate('wildcard','*RECORD*')) == 1 % search for RECORD
                 sCount = sCount+1; % document number of starts
                 OsRow(sCount) = iD; % start row for original data array at RECORD for each start
             end
        elseif strcmp(Odata{iD,1},'END') && ~isempty(Odata{iD,1}) == 1
            eCount = eCount+1; % document number of ends
            OeRow(eCount) = iD; % end row for original data array at END for each end
        end
    end
end

% determine number of blocks ofg trials
if sCount>=eCount  % are there at least as many starts as ends?
    nBlocks = eCount; % set the number of blocks of trails in the file
elseif sCount<eCount % are there fewer starts than ends?
    nBlocks = sCount; % set the number of blocks of trials in the file
end

data = cell(length(Odata),9); % preallocate size of data
data = repmat({data},1,nBlocks);

for B = 1:nBlocks; % cycle through blocks of trials
    % create block structure
    data{1,B} = Odata(OsRow(B):OeRow(B),1:9); % assign only relevant rows to data array
    % preallocate block structure fields and sizes
    timepts = cell(length(data{B}),4);
    timepts(1,1:4) = {'time', 'X coordinate', 'Y coordinate', 'pupil size'};
    timepts = repmat({timepts},1,nBlocks);
    events = cell(length(data{B}),2);
    events(1,1:2) = {'time', 'event name'};
    events = repmat({events},1,nBlocks);
    blink = cell(length(data{B}),3); 
    blink(1,1:3) = {'start time', 'end time', 'duration'};
    blink = repmat({blink},1,nBlocks);
    fix = cell(length(data{B}),6);
    fix(1,1:6) = {'start time', 'end time', 'angular resolution', 'average X position', 'average Y position', 'average pupil size'};
    fix = repmat({fix},1,nBlocks);
    sacc = cell(length(data{B}),9);
    sacc(1,1:9) = {'start time', 'end time', 'duration', 'start X position', 'start Y position', 'end X position', 'end Y position', 'amplitude (degrees)', 'velocity (degrees/second)'};
    sacc = repmat({sacc},1,nBlocks);
end
    
block = struct('timepts',timepts,'events',events,'blink',blink,'fix',fix,'sacc',sacc); % set output structure fields
                
%% Assign events and timepoints

tRow(1:nBlocks) = 2; %set up rows for indexing into timepts
eventRow(1:nBlocks) = 2; % set up event count variable 
blinkRow(1:nBlocks) = 2; % set up blink count variable 
fixRow(1:nBlocks) = 2; % set up fixation count variable
saccRow(1:nBlocks) = 2; % set up saccade count variable

for iB = 1:nBlocks; % cycle through blocks of trials
        % find positions in block iB of data array where events occur
        clear MSGarray posEvents
        MSGarray = cell(size(data{iB}(1:end,1))); % create array of str MSG the length of data{iB}
        MSGarray(cellfun('isempty',MSGarray)) = {'MSG'};
        posEvents = find(cellfun(@strcmp,data{iB}(:,1),MSGarray)==1); % find the position of MSG string thoughout first column of data{iB}
        if ~isempty(posEvents)
            for iE = 1:length(posEvents); % cycle through events
                clear S is_number iS
                S = data{iB}{posEvents(iE),2}; % set S as string to examine
                is_number = zeros(1,length(S));
                for iS = 1:length(S); % cycle through positions within string
                    is_number(iS) = isnumeric(str2double(S(iS))) && ~isnan(str2double(S(iS))) && str2double(S(iS))~=1i; % create index of (non 1i) numeric characters
                end
                posEventNmStart = find(is_number==0,1); % find the position of the first non-number in S
                if ~isempty(posEventNmStart)
                    block(iB).events(eventRow(iB),2) = {S(posEventNmStart+1:end)}; % set event label excluding preceding time stamp
                    block(iB).events(eventRow(iB),1) = {S(1:posEventNmStart-1)}; % set event time
                    eventRow(iB) = eventRow(iB)+1;
                end
            end
        end
    for iL = 1:length(data{iB}); % cycle through lines of data cell array
        if isnumeric(str2double(data{iB}{iL,1})) && ~isnan(str2double(data{iB}{iL,1})) == 1 % search for numbers in the first column of data
            block(iB).timepts(tRow(iB),1) = (data{iB}(iL,1)); % store timepoint
            block(iB).timepts(tRow(iB),2) = (data{iB}(iL,2)); % store X coordinate
            block(iB).timepts(tRow(iB),3) = (data{iB}(iL,3)); % store Y coordinate
            block(iB).timepts(tRow(iB),4) = (data{iB}(iL,4)); % store pupil size
            tRow(iB) = tRow(iB)+1;
%         elseif any(posEvents==iL) % is iL the location of an event MSG?
%             clear S is_number iS
%             S = data{iB}{iL,2}; % set S as string to examine
%             is_number = zeros(1,length(S));
%             for iS = 1:length(S); % cycle through positions within string
%                 is_number(iS) = isnumeric(str2double(S(iS))) && ~isnan(str2double(S(iS))) && str2double(S(iS))~=1i || isspace(S(iS)); % create index of (non 1i) numeric characters and spaces
%             end
%             posEventNmStart = find(is_number==0,1); % find the position of the first non-number in S
%             if ~isempty(posEventNmStart)
%                 block(iB).timepts(tRow(iB)+1:end,6) = {S(posEventNmStart:end)}; % set event label excluding preceding time stamp
%                 disp(['documented event: ' S(posEventNmStart:end) ' at tRow: ' num2str(tRow(iB)+1) ' and data row: ' num2str(iL)]) % test
%                 block(iB).timepts = block(iB).timepts([1:tRow(iB)-1 tRow(iB)+1:end],1:6); % remove extra row, which marks event in data, from timepts array
%             end
        elseif regexp(data{iB}{iL,1},regexptranslate('wildcard','SBLINK R*')) == 1; % search for blinks
            block(iB).blink(blinkRow(iB),1) = data{iB}(iL+1,1); % document time point of blink start
%             block(iB).timepts(tRow(iB)+1:end,5) = {nan}; % document blink
%             block(iB).timepts = block(iB).timepts([1:tRow(iB)-1 tRow(iB)+1:end],1:6); % remove extra row, which marks event in data, from timepts array
        elseif regexp(data{iB}{iL,1},regexptranslate('wildcard','SFIX R*')) == 1 % search for fixations
            block(iB).fix(fixRow(iB),1) = data{iB}(iL+1,1); % document time point of fixation start
%             block(iB).timepts(tRow(iB)+1:end,5) = {1}; % document fixation
%             block(iB).timepts = block(iB).timepts([1:tRow(iB)-1 tRow(iB)+1:end],1:6); % remove extra row, which marks event in data, from timepts array
        elseif regexp(data{iB}{iL,1},regexptranslate('wildcard','SSACC R*')) == 1 % search for saccades
            block(iB).sacc(saccRow(iB),1) = data{iB}(iL+1,1); % document time point of saccade start
%             block(iB).timepts(tRow(iB)+1:end,5) = {0}; % document saccade
%             block(iB).timepts = block(iB).timepts([1:tRow(iB)-1 tRow(iB)+1:end],1:6); % remove extra row, which marks event in data, from timepts array
        elseif regexp(data{iB}{iL,1},regexptranslate('wildcard','EBLINK R*')) == 1 % search for blinks ending
            block(iB).blink(blinkRow(iB),2:3) = data{iB}(iL,2:3); % document blink end time point & duration
            blinkRow(iB) = blinkRow(iB)+1; % add one to blink count
%             block(iB).timepts = block(iB).timepts([1:tRow(iB)-1 tRow(iB)+1:end],1:6); % remove extra row, which marks event in data, from timepts array
        elseif regexp(data{iB}{iL,1},regexptranslate('wildcard','EFIX R*')) == 1 % search for fixations ending
            block(iB).fix(fixRow(iB),2:6) = data{iB}(iL,2:6); % document fixation end time point, angular resolution, avg X, avg Y, & avg pupil
            fixRow(iB) = fixRow(iB)+1; % add one to blink count
%             block(iB).timepts = block(iB).timepts([1:tRow(iB)-1 tRow(iB)+1:end],1:6); % remove extra row, which marks event in data, from timepts array
        elseif regexp(data{iB}{iL,1},regexptranslate('wildcard','ESACC R*')) == 1 % search for saccades ending
            block(iB).sacc(saccRow(iB),2:9) = data{iB}(iL,2:9); % document saccade end time point, duration, start X, start Y, end X, end Y, amp, velocity
            saccRow(iB) = saccRow(iB)+1; % add one to blink count
%             block(iB).timepts = block(iB).timepts([1:tRow(iB)-1 tRow(iB)+1:end],1:6); % remove extra row, which marks event in data, from timepts array
        end
    end

    % resize funtion output
    tRow(iB) = tRow(iB)-1;
    eventRow(iB) = eventRow(iB)-1;
    blinkRow(iB) = blinkRow(iB)-1;
    fixRow(iB) = fixRow(iB)-1;
    saccRow(iB) = saccRow(iB)-1;

    block(iB).timepts = block(iB).timepts(1:tRow(iB),1:4);
    block(iB).events = block(iB).events(1:eventRow(iB),1:2);
    block(iB).blink = block(iB).blink(1:blinkRow(iB),1:3);
    block(iB).fix = block(iB).fix(1:fixRow(iB),1:6);
    block(iB).sacc = block(iB).sacc(1:saccRow(iB),1:9);
end

end

