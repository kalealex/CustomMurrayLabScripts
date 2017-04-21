function [accuracy,resptimes,timepts,pulseCount] = AK_GABAinASD_logfileParse(filename,startRow,endRow)
%AK_GABAinASD_logfileParse takes a log file name as input and returns
%response accuracy, response times, a cell array of timepoints, and a pulse count for the logfile.
%   imports trial#s, event types, codes, and times from tab-delimited log file
%   parses 'greencirc', 'other' non-stimulus shapes, and 'Response' events at timepoints
%   calculates response times; returns structure with vectors including and
%   excluing nonhit responstimes
%   determines hits, misses, and errors
%   determines response accuracy as hit rate, miss rate, and false alarm
%   rate overall (includes # presentations of green circles, other shapes
%   and # of responses in accuracy output)
%   returns list of response times for hits only and all response times in
%   response structure
%   assigns conditions and accuracy information by stimulus and response
%   events in timepts array


if ischar(filename) == 0; % check that filename is a string
    disp('filename is not a string')
end

% import events and times from tab-delimited log file
%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 4;
    endRow = inf;
end

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


% parse data from cell array
%% create Cond structure to identify condition per trial based on logfile name

% designate identifying wildcard strings for each potential type of logfile
MTloc_str = '*MTlocalizer*';
Contrast_str = '*contrast*';
Supp_str = '*suppression*';
Summ_str = '*summation*';
V1_str = '*V1localizer*';

% preallocate fields of structure
cond_name = cell(1);
bmp_name = cell(1);
trial_numbers = nan(500,1); % set nans = [] at end of function
% hits = 0;
% misses = 0;
% false_alarms = 0;
% response_times = cell(1);
% stim_presentations = 0;
% other_presentations = 0;

% generate cond struct w/ above fields
if AK_findStrMatch({filename},MTloc_str) == 1;
    Cond(1:2) = struct('label',cond_name,'bmp',bmp_name,'trials',trial_numbers);
    Cond(1).label = 'moving';
    Cond(1).bmp = '*loc5_orient*.bmp'; % identifying bmp is the same for moving and static; moving is new orientation each trial
    Cond(2).label = 'static';
    Cond(2).bmp = '*loc5_orient*.bmp';
elseif AK_findStrMatch({filename},Contrast_str) == 1;
    Cond(1:3) = struct('label',cond_name,'bmp',bmp_name,'trials',trial_numbers);
    Cond(1).label = 'fix';
    Cond(1).bmp = 'blank*.bmp';
    Cond(2).label = 'hi';
    Cond(2).bmp = '*size1_contrast2*.bmp';
    Cond(3).label = 'low';
    Cond(3).bmp = '*size1_contrast1*.bmp';
elseif AK_findStrMatch({filename},Supp_str) == 1;
    Cond(1:2) = struct('label',cond_name,'bmp',bmp_name,'trials',trial_numbers);
    Cond(1).label = 'big';
    Cond(1).bmp = '*size2*.bmp';
    Cond(2).label = 'small';
    Cond(2).bmp = '*size1*.bmp';
elseif AK_findStrMatch({filename},Summ_str) == 1;
    Cond(1:2) = struct('label',cond_name,'bmp',bmp_name,'trials',trial_numbers);
    Cond(1).label = 'big';
    Cond(1).bmp = '*size2*.bmp';
    Cond(2).label = 'small';
    Cond(2).bmp = '*size1*.bmp';
elseif AK_findStrMatch({filename},V1_str) == 1;
    Cond(1:2) = struct('label',cond_name,'bmp',bmp_name,'trials',trial_numbers);
    Cond(1).label = 'surr';
    Cond(1).bmp = '*check*';
    Cond(2).label = 'tar'; 
    Cond(2).bmp = '*check*';
end

%% parse 'greencirc' and 'Response' at timepoints and assign trial numbers to conditions

L = length(data(:,1)); % length of imported cell array
presentations_stim = 0; % create variable for number of stimulus presentations
presentations_other = 0; % create variable for number of non-stimulus shape presentations
responses = 0; % create variable for number of participant responses
timepts = cell(L,6); % set up timepoints cell array
timepts(1,1:6) = {'trial','event','time','resptime','accuracy','condition'}; % timepoints cell array labels
row = 2; % set up row variable for indexing into timepoints cell array
condTrial_count = zeros(1,length(Cond)); % set up counter for indexing into trials for each condition
last_orient = []; % orientation comparison string for MT localizers only
pulseCount = 0; % create pulse counter

for i = 1:L; % run through cell array
    
    if strcmp(Cond(1).label,'moving') == 1; % is the logfile an MT localizer?
        if AK_findStrMatch(data(i,3),Cond(1).bmp) == 1; % does the code match the bmp name? (both conditions have the same Cond.bmp)
            if isempty(last_orient) == 1; % is it the first trial?
                last_orient = str2double(data{i,3}(20)); % document last orientation
                this_orient = str2double(data{i,3}(20)); % document current orientation
            else
                clear this_orient
                this_orient = str2double(data{i,3}(20)); % document current orientation
            end
            if last_orient == this_orient % is this orientation the same as the last orientation?
                condTrial_count(2) = condTrial_count(2)+1; % add to trial count index for static cond
                Cond(2).trials(condTrial_count(2)) = data{i,1}; % document trial number
            else % if this bmp is different than the last one
                condTrial_count(1) = condTrial_count(1)+1; % add to trial count index for moving cond
                Cond(1).trials(condTrial_count(1)) = data{i,1}; % document trial number
            end
            clear last_orient
            last_orient = this_orient; % set last orientation to this orientation
        end 
    elseif strcmp(Cond(1).label,'surr') == 1; % is the logfile a V1 localizer?
        if AK_findStrMatch(data(i,3),Cond(1).bmp) == 1; % does the code match the bmp name? (both conditions have the same Cond.bmp)
            if length(data{i,3}) > 5 && strcmp(data{i,3}(6),'s') == 1; % is the 6th letter of the string an 's'? Is the condition named 'checks*' or 'check*'?
                condTrial_count(1) = condTrial_count(1)+1; % add to trial count index for surround cond
                Cond(1).trials(condTrial_count(1)) = data{i,1}; % document trial number
            else % if the string does not have an 's' in the 6th position
                condTrial_count(2) = condTrial_count(2)+1; % add to trial count index for target cond
                Cond(2).trials(condTrial_count(2)) = data{i,1}; % document trial number
            end
        end
    else % if the logfile is not an MT or V1 localizer
        for c = 1:length(Cond); % cycle through conditions for this logfile
            if AK_findStrMatch(data(i,3),Cond(c).bmp) == 1; % does the code match the bmp name?
                condTrial_count(c) = condTrial_count(c)+1; % add to trial count index for cond c
                Cond(c).trials(condTrial_count(c)) = data{i,1}; % document trial number
            end
        end
    end 
        
    if strcmp(data(i,3),'greencirc') == 1; % search for greencirc
        timepts(row,1) = data(i,1); % document trial number
        timepts(row,2) = {'greencirc'}; % document event
        presentations_stim = presentations_stim +1; % count number of stimulus presentations
        timepts(row,3) = data(i,4); % document timepoint
        row = row+1; % move to next row of timepoints
    elseif strcmp(data(i,3),'other') == 1; % search for greencirc
        timepts(row,1) = data(i,1); % document trial number
        timepts(row,2) = {'other'}; % document event
        presentations_other = presentations_other +1; % count number of non-stimulus shape presentations
        timepts(row,3) = data(i,4); % document timepoint
        row = row+1; % move to next row of timepoints
    elseif strcmp(data(i,2),'Response') == 1; % search for Response
        timepts(row,1) = data(i,1); % document trial number
        timepts(row,2) = {'Response'}; % document event
        responses = responses+1; % count number of participant responses
        timepts(row,3) = data(i,4); % document timepoint
        row = row+1; % move to next row of timepoints
    elseif strcmp(data(i,2),'Pulse') == 1; % search for Pulse
        pulseCount = pulseCount+1;
    end
end

row = row-1; % adjust row count

%% determine response times, accuracy, and condition labels

hits = 0; % set up variable to count hits
misses = 0; % set up variable to count misses
false_alarms = 0; % set up variable to count errors

for i = 1:row; % run through timepoints
    
    for c = 1:length(Cond); % cycle through conditions for this logfile
        if ~isempty(intersect(timepts{i,1},Cond(c).trials)) && i ~= 1 % is the trial for row i of timepts part of the trial list for condition c?
            timepts(i,6) = {Cond(c).label};
        end
    end
    
    if strcmp(timepts(i,2),'greencirc') == 1; % is the event a stim presentation?
        if strcmp(timepts(i+1,2),'Response') == 1; % is the stim presentation followed by a response?
            timepts(i+1,4) = {cell2mat(timepts(i+1,3))-cell2mat(timepts(i,3))}; % document response time
            if cell2mat(timepts(i+1,3))-cell2mat(timepts(i,3)) < 20000 % is response time less than 2 seconds?
                timepts(i,5) = {'hit'}; % yes: document hit
                hits = hits +1; % count the number of hits
            else timepts(i,5) = {'miss'}; % no: document miss
                timepts(i+1,5) = {'nonhit'}; % mark response time as nonhit for indexing (resptimes.hits)
                misses = misses+1; % count the number of misses
            end
        else timepts(i,5) = {'miss'}; % no: document miss
            misses = misses+1; % count the number of misses
        end
    elseif strcmp(timepts(i,2),'other') == 1; % is the event a non-stim presentation? 
        if strcmp(timepts(i+1,2),'Response') == 1; % is the non-stim presentation followed by a response?
            timepts(i+1,4) = {cell2mat(timepts(i+1,3))-cell2mat(timepts(i,3))}; % document response time
            timepts(i+1,5) = {'false alarm'}; % document false alarm response
            false_alarms = false_alarms+1; % count the number of false alarms
        end
    elseif strcmp(timepts(i,2),'Response') == 1; % is the event a response?
%         timepts(i,1) = timepts(i-1,1); % change trial number of response so that the response is associated with the preceding event (in order to separate response times by conditions)
        if strcmp(timepts(i-1,2),'greencirc') == 0; % is the response NOT preceded by a stim?
            timepts(i,5) = {'false alarm'}; % document false alarm response
            false_alarms = false_alarms+1; % count the number of false alarms
        end
    end
end

timepts = timepts(1:row,:); % remove empty cells at the end of timepts array

%% response times structure containing vectors for all response times as well as response times for hits only

resptimes.all = timepts(2:end,4); % extract response times from timepoints cell array

notNonhitIndex = find(~strcmp(timepts(2:end,5),'nonhit')); % create index for responses which are not nothits
notFalseAlarmIndex = find(~strcmp(timepts(2:end,5),'false alarm')); % create index for responses which are not errors
hitIndex = intersect(notNonhitIndex,notFalseAlarmIndex); % create index of responses which are hits (not errors or nonhits)
resptimes.hits = resptimes.all(hitIndex); % filter out response times for nonhits and errors

resptimes.all = resptimes.all(~cellfun('isempty',resptimes.all)); % remove empty cells
resptimes.hits = resptimes.hits(~cellfun('isempty',resptimes.hits));

resptimes.all = cell2mat(resptimes.all); % transform cell arrays into structures of vectors
resptimes.hits = cell2mat(resptimes.hits);

if isempty(resptimes.all); % account for files with no responses
    resptimes.all = nan;
end

if isempty(resptimes.hits); % account for files with no responses
    resptimes.hits = nan;
end

%% accuracy metrics

hitRate = hits/presentations_stim; % calculate hit rate
missRate = misses/presentations_stim; %calculate miss rate
false_alarmRate = false_alarms/presentations_other; % calculate error rate
accuracy = [hitRate missRate false_alarmRate responses presentations_stim presentations_other]; % create vector of accuracy metrics to summarize participant performance

end

