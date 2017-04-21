function [accuracy,resptimes,timepts,pulseCount] = AK_SupSum_logfileParse(file_dir)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ischar(file_dir) == 0; % check that filename is a string
    disp('filename is not a string')
end
if nargin<1
    file_dir = input('What full file directory (string) would you like to analyze?');
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
fileID = fopen(file_dir,'r');

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

%% Create output variable

dataArray([1, 4]) = cellfun(@(x) num2cell(x), dataArray([1, 4]), 'UniformOutput', false);
data = [dataArray{1:end-1}];

%% parse data from cell array

% setup
L = length(data(:,1)); % length of imported cell array
presentations_stim = 0; % create variable for number of stimulus presentations
presentations_other = 0; % create variable for number of noise color presentations
responses = 0; % create variable for number of participant responses
pulseCount = 0; % create pulse counter
timepts = cell(L,5); % set up timepoints cell array
timepts(1,1:5) = {'trial','event','time','resptime','performance'}; % timepoints cell array labels
row = 2; % set up row variable for indexing into timepoints cell array

% run through rows of data array filing in timepts: trial, event & time
% fields; document number of signal and noise trials, responses and number of pulses
for iD = 1:L;
    if strcmp(data(iD,3),'greencirc') == 1; % search for greencirc
        timepts(row,1) = data(iD,1); % document trial number
        timepts(row,2) = {'greencirc'}; % document event
        presentations_stim = presentations_stim +1; % count number of stimulus presentations
        timepts(row,3) = data(iD,4); % document timepoint
        row = row+1; % move to next row of timepoints
    elseif strcmp(data(iD,3),'other') == 1; % search for noise presentations
        timepts(row,1) = data(iD,1); % document trial number
        timepts(row,2) = {'other'}; % document event
        presentations_other = presentations_other +1; % count number of noise color presentations
        timepts(row,3) = data(iD,4); % document timepoint
        row = row+1; % move to next row of timepoints
    elseif strcmp(data(iD,2),'Response') == 1; % search for Response
        timepts(row,1) = data(iD,1); % document trial number
        timepts(row,2) = {'Response'}; % document event
        responses = responses+1; % count number of participant responses
        timepts(row,3) = data(iD,4); % document timepoint
        row = row+1; % move to next row of timepoints
    elseif strcmp(data(iD,2),'Pulse') == 1; % search for Pulse
        pulseCount = pulseCount+1;
    end
end

% adjust row count
row = row-1; 

% setup
hits = 0; % set up variable to count hits
misses = 0; % set up variable to count misses
false_alarms = 0; % set up variable to count errors

% fill in timepts: resptimes & performance fields; count hits misses and
% false alarms
for iT = 1:row
    if strcmp(timepts(iT,2),'greencirc') == 1; % is the event a stim presentation?
        if strcmp(timepts(iT+1,2),'Response') == 1; % is the stim presentation followed by a response?
            timepts(iT+1,4) = {cell2mat(timepts(iT+1,3))-cell2mat(timepts(iT,3))}; % document response time
            if cell2mat(timepts(iT+1,3))-cell2mat(timepts(iT,3)) < 20000 % is response time less than 2 seconds?
                timepts(iT,5) = {'hitStim'}; % yes: document hit
                timepts(iT+1,5) = {'hitResp'}; % yes: document hit
                hits = hits +1; % count the number of hits
            else
                timepts(iT,5) = {'missStim'}; % no: document miss
                timepts(iT+1,5) = {'nonhitResp'}; % mark response time as nonhit for indexing (resptimes.hits)
                misses = misses+1; % count the number of misses
            end
        else
            timepts(iT,5) = {'miss'}; % no: document miss
            misses = misses+1; % count the number of misses
        end
    elseif strcmp(timepts(iT,2),'other') == 1; % is the event a noise presentation? 
        if strcmp(timepts(iT+1,2),'Response') == 1; % is the noise presentation followed by a response?
            timepts(iT+1,4) = {cell2mat(timepts(iT+1,3))-cell2mat(timepts(iT,3))}; % document response time
            timepts(iT+1,5) = {'falseAlarmResp'}; % document false alarm response
            false_alarms = false_alarms+1; % count the number of false alarms
        end
    elseif strcmp(timepts(iT,2),'Response') == 1; % is the event a response?
        if strcmp(timepts(iT-1,2),'greencirc') == 0; % is the response NOT preceded by a stim?
            timepts(iT,5) = {'falseAlarmResp'}; % document false alarm response
            false_alarms = false_alarms+1; % count the number of false alarms
        end
    end
end

% remove empty cells at the end of timepts array
timepts = timepts(1:row,:); 

% create resptimes output:
% separate hit, nonhit & false alarm responses with indices for each
hitRespIdx = strcmp(timepts(:,5),'hitResp');
nonhitRespIdx = strcmp(timepts(:,5),'nonhitResp');
falseAlarmRespIdx = strcmp(timepts(:,5),'falseAlarmResp');
% extract response times from timepoints cell array
resptimes.all = timepts(2:end,4);
resptimes.hits = timepts(hitRespIdx,4);
resptimes.nonhits = timepts(nonhitRespIdx,4);
resptimes.falseAlarms = timepts(falseAlarmRespIdx,4);
% remove empty cells
resptimes.all = resptimes.all(~cellfun(@isempty,resptimes.all));
resptimes.hits = resptimes.hits(~cellfun(@isempty,resptimes.hits));
resptimes.nonhits = resptimes.nonhits(~cellfun(@isempty,resptimes.nonhits));
resptimes.falseAlarms = resptimes.falseAlarms(~cellfun(@isempty,resptimes.falseAlarms));
% convert to matrices
resptimes.all = cell2mat(resptimes.all)./10000;
resptimes.hits = cell2mat(resptimes.hits)./10000;
resptimes.nonhits = cell2mat(resptimes.nonhits)./10000;
resptimes.falseAlarms = cell2mat(resptimes.falseAlarms)./10000;
% account for files with no responses
if isempty(resptimes.all); 
    resptimes.all = nan;
end
if isempty(resptimes.hits); 
    resptimes.hits = nan;
end
if isempty(resptimes.nonhits); 
    resptimes.nonhits = nan;
end
if isempty(resptimes.falseAlarms); 
    resptimes.falseAlarms = nan;
end

% create accuracy output:
accuracy.hitRate = hits/presentations_stim; % calculate hit rate
accuracy.missRate = misses/presentations_stim; %calculate miss rate
accuracy.falseAlarmRate = false_alarms/presentations_other; % calculate false alarm rate
accuracy.nResponses = responses;
accuracy.nSignalTrials = presentations_stim; 
accuracy.nNoiseTrials = presentations_other;


end

