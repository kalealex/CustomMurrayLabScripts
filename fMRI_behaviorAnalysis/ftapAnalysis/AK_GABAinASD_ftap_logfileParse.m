function [accuracy,resptimes,timepts] = AK_GABAinASD_ftap_logfileParse(filename,startRow,endRow)
%AK_GABAinASD_logfileParse takes a log file name as input and returns
%response accuracy, response times, and a cell array of timepoints.
%   imports trial#s, event types, codes, and times from tab-delimited log file
%   parses 'standard', 'variable', and 'Response' events at timepoints
%   calculates response times; returns structure with vectors including and
%   excluing nothit responstimes, separated by type of stimulus event
%   determines hits and misses per type of stimulus event, as well as errors
%   determines response accuracy as hit rate and miss rate per type of stimulus event, as well as error rate


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
%% parse 'greencirc' and 'Response' at timepoints

L = length(data(:,1)); % length of imported cell array
presentationsStandard = 0; % create variable for number of stimulus presentations
presentationsVariable = 0; % create variable for number of stimulus presentations
responses = 0; % create variable for total number of participant responses, regardless of stimulus
timepts = cell(L,5); % set up timepoints cell array
timepts(1,1:5) = {'trial','event','time','resptime','accuracy'}; % timepoints cell array labels
row = 2; % set up row variable for indexing into timepoints cell array

for i = 1:L; % run through cell array
    if strcmp(data(i,3),'standard') == 1; % search for standard
        timepts(row,1) = data(i,1); % document trial number
        timepts(row,2) = {'standard'}; % document event
        presentationsStandard = presentationsStandard +1; % count number of standard stimulus presentations
        timepts(row,3) = data(i,4); % document timepoint
        row = row+1; % move to next row of timepoints
    elseif strcmp(data(i,3),'variable') == 1; % search for variable
        timepts(row,1) = data(i,1); % document trial number
        timepts(row,2) = {'variable'}; % document event
        presentationsVariable = presentationsVariable +1; % count number of variable stimulus presentations
        timepts(row,3) = data(i,4); % document timepoint
        row = row+1; % move to next row of timepoints
    elseif strcmp(data(i,2),'Response') == 1; % search for Response
        timepts(row,1) = data(i,1); % document trial number
        timepts(row,2) = {'Response'}; % document event
        responses = responses+1; %count total number of participant responses
        timepts(row,3) = data(i,4); % document timepoint
        row = row+1; % move to next row of timepoints
    end
end

row = row-1; % adjust row count

%% determine response times and accuracy

hitsStandard = 0; % set up variable to count hits for standard stim
hitsVariable = 0; % set up variable to count hits for variable stim
missesStandard = 0; % set up variable to count misses for standard stim
missesVariable = 0; % set up variable to count misses for variable stim
errors = 0; % set up variable to count errors for stim

for i = 1:row; % run through timepoints
    if strcmp(timepts(i,2),'standard') == 1; % is the event a standard stim presentation?
        if strcmp(timepts(i+1,2),'Response') == 1; % is the stim presentation followed by a response?
            timepts(i+1,4) = {cell2mat(timepts(i+1,3))-cell2mat(timepts(i,3))}; % document response time
            if cell2mat(timepts(i+1,3))-cell2mat(timepts(i,3)) < 10000 % is response time less than 1 second?
                timepts(i,5) = {'hit'}; % yes: document hit
                timepts(i+1,5) = {'hitStd'}; % mark response as hit for standard stim
                hitsStandard = hitsStandard +1; % count the number of hits for standard stim
            else timepts(i,5) = {'miss'}; % no: document miss
                timepts(i+1,5) = {'nonhitStd'}; % mark response time as nonhit for standard stim
                missesStandard = missesStandard+1; % count the number of misses for standard stim
            end
        else timepts(i,5) = {'miss'}; % no: document miss
            missesStandard = missesStandard+1; % count the number of misses for standard stim
        end
    elseif strcmp(timepts(i,2),'variable') == 1; % is the event a variable stim presentation?
        if strcmp(timepts(i+1,2),'Response') == 1; % is the stim presentation followed by a response?
            timepts(i+1,4) = {cell2mat(timepts(i+1,3))-cell2mat(timepts(i,3))}; % document response time
            if cell2mat(timepts(i+1,3))-cell2mat(timepts(i,3)) < 10000 % is response time less than 1 second?
                timepts(i,5) = {'hit'}; % yes: document hit
                timepts(i+1,5) = {'hitVar'}; % mark response as hit for variable stim
                hitsVariable = hitsVariable +1; % count the number of hits for variable stim
            else timepts(i,5) = {'miss'}; % no: document miss
                timepts(i+1,5) = {'nonhitVar'}; % mark response time as nonhit for variable stim
                missesVariable = missesVariable+1; % count the number of misses for variable stim
            end
        else timepts(i,5) = {'miss'}; % no: document miss
            missesVariable = missesVariable+1; % count the number of misses for varibale stim
        end    
    elseif strcmp(timepts(i,2),'Response') == 1; % is the event a response?
        if strcmp(timepts(i-1,2),'standard') == 0 && strcmp(timepts(i-1,2),'variable') == 0; % is the response NOT preceded by a stim?
            timepts(i,5) = {'error'}; % document erroneous response
            errors = errors+1; % count the number of errors
        end
    end
end

timepts = timepts(1:row,:); % remove empty cells at the end of timepts array

%% response times structure containing vectors for all response times as well as response times for hits only

resptimes.all = timepts(2:end,4); % extract response times from timepoints cell array

StdIndex = ~cellfun(@isempty,regexp(timepts(2:end,5),regexptranslate('wildcard','*Std'))); % create index for all responses to standard stim
VarIndex = ~cellfun(@isempty,regexp(timepts(2:end,5),regexptranslate('wildcard','*Var'))); % create index for all responses to variable stim
hitStdIndex = strcmp(timepts(2:end,5),'hitStd'); % create index for responses which are hits for standard stim
hitVarIndex = strcmp(timepts(2:end,5),'hitVar'); % create index for responses which are hits for varibale stim

resptimes.allStd = resptimes.all(StdIndex); % response times for standard stim
resptimes.allVar = resptimes.all(VarIndex); % response times for variable stim
resptimes.hitsStd = resptimes.all(hitStdIndex); % response times which are hits for standard stim
resptimes.hitsVar = resptimes.all(hitVarIndex); % response times which are hits for variable stim

resptimes.all = resptimes.all(~cellfun('isempty',resptimes.all)); % remove empty cells
resptimes.allStd = resptimes.allStd(~cellfun('isempty',resptimes.allStd));
resptimes.allVar = resptimes.allVar(~cellfun('isempty',resptimes.allVar));
resptimes.hitsStd = resptimes.hitsStd(~cellfun('isempty',resptimes.hitsStd));
resptimes.hitsVar = resptimes.hitsVar(~cellfun('isempty',resptimes.hitsVar));

resptimes.all = cell2mat(resptimes.all); % transform cell arrays into structures of vectors
resptimes.allStd = cell2mat(resptimes.allStd);
resptimes.allVar = cell2mat(resptimes.allVar);
resptimes.hitsStd = cell2mat(resptimes.hitsStd);
resptimes.hitsVar = cell2mat(resptimes.hitsVar);

if isempty(resptimes.all); % account for files with no responses
    resptimes.all = nan;
end

if isempty(resptimes.allStd); % account for files with no responses
    resptimes.allStd = nan;
end

if isempty(resptimes.allVar); % account for files with no responses
    resptimes.allVar = nan;
end

if isempty(resptimes.hitsStd); % account for files with no responses
    resptimes.hitsStd = nan;
end

if isempty(resptimes.hitsVar); % account for files with no responses
    resptimes.hitsVar = nan;
end

%% accuracy metrics

hitRateStd = hitsStandard/presentationsStandard; % calculate hit rate for standard stim
hitRateVar = hitsVariable/presentationsVariable; % calculate hit rate for variable stim
missRateStd = missesStandard/presentationsStandard; %calculate miss rate for standard stim
missRateVar = missesVariable/presentationsVariable; %calculate miss rate for variable stim
errorRate = errors/responses; % calculate error rate
accuracy = [hitRateStd hitRateVar missRateStd missRateVar errorRate responses]; % create vector of accuracy metrics to summarize participant performance

end

