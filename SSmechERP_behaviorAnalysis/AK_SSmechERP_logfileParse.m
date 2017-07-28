function [ overall, stereoAlign, C1, SSmech ] = AK_SSmechERP_logfileParse( subj_initials, data_dir )
%AK_SSmechERP_logfileParse parses the .log files containing behavioral data
%from the SSmech ERP experiment run by Michael-Paul Schallmo and Alex Kale
%in the Murray Lab (UW Psychology) in the summer of 2016. Given a
%participant's initials and a data directory, the function will retrieve
%and analyze the behavioral data for the participant, separating data based
%on whether it is from the stereoscope alignment task, the C1 test one-back
%task, or the SSmech plaid detection task.
%   INPUT:
%       subj_initials: a subject's initials (as a string) exactly as they
%           appear in the .log file names; if a subject's initials do not
%           match on all files, they must be manually changed so that they
%           match in order for this function to group the data properly
%       data_dir: the directery where the data is saved (as a sting)
%   OUTPUT:
%       overall: a data structure containing summary statistics about this
%           subject's overall performance in the experiment; fields are as
%           follows:
%           finalAccuracyStereoAlign: lists the accuracy the participant
%               attained on the final attempt of the stereoscope alignment
%               task for each file of this type
%           avgAccuracyC1: the average accuracy the participant attained
%               across sessions of the C1 test one-back task
%           avgHitRateSSmech: the average hit rate the participant attained
%               across sessions of the SSmech plaid detection task
%           avgMissRateSSmech: the average miss rate the participant
%               attained across sessions of the SSmech plaid detection task
%           avgFalseAlarmRateSSmech: the average false alarm rate the
%               participant attained across sessions of the SSmech plaid
%               detection task 
%           avgCorrectRejectionRateSSmech: the average correct rejection
%               rate the participant attained across sessions of the SSmech
%               plaid detection task
%       stereoAlign: a data structure containing accuracies on the
%           stereoscope alignment task for each file of this type among the
%           subject's data; there is one unit in the structure array per
%           file; fields are as follows: 
%           filename: the .log file name for each file of this type
%           date: the date and time that each file was written
%           accuracy: the accuracy attained by the participant on each
%               attempt at the stereoscope alignment task for each file
%       C1: a data structure containing accuracies on the C1 test one-back
%           task for each file of this type among the subject's data; there
%           is one unit in the structure array per file; fields are as
%           follows:
%           filename: the .log file name for each file of this type
%           date: the date and time that each file was written
%           correct: a cell array where each cell corresponds to a block
%               and contains a logical array with correct (1) and incorrect
%               (0) scoring for each trial
%           noResponse: a cell array where each cell corresponds to a block
%               and contains a logical array with no response (1) and
%               response (0) recorded for each trial
%           accuracy: an array of accuracies attained in each block of the
%               C1 test one-back task 
%           accuracyAcrossBlocks: the overall accuracy across blocks within
%               each file
%       SSmech: a data structure containing hit, miss, false alarm, and
%           correct rejection data on the SSmech plaid detection task for
%           each file of this type among the subject's data; there 
%           is one unit in the structure array per file; fields are as
%           follows:
%           filename: the .log file name for each file of this type
%           date: the date and time that each file was written
%           hits: a cell array where each cell corresponds to a block
%               and contains a logical array with hit (1) and miss (0)
%               scoring for each trial
%           misses: a cell array where each cell corresponds to a block
%               and contains a logical array with miss (1) and hit (0)
%               scoring for each trial
%           falseAlarms: a cell array where each cell corresponds to a
%               block and contains a logical array with false alarm (1) and
%               correct rejection (0) scoring for each trial
%           correctRejections: a cell array where each cell corresponds to a
%               block and contains a logical array with correct rejection
%               (1) and false alarm (0) scoring for each trial
%           isSignal: a cell array where each cell corresponds to a
%               block and contains a logical array with plaid (1) and
%               non-plaid (0) stimuli documented for each trial
%           hitRate: an array of hit rates attained in each block of the
%               SSmech plaid detections task
%           missRate: an array of miss rates attained in each block of the
%               SSmech plaid detections task
%           falseAlarmRate: an array of false alarm rates attained in each
%           block of the SSmech plaid detections task
%           correctRejectionRate: an array of correct rejection rates
%           attained in each block of the SSmech plaid detections task
%           hitRateAcrossBlocks: the average hit rate across blocks of the
%               experiment within each file
%           missRateAcrossBlocks: the average miss rate across blocks of the
%               experiment within each file
%           falseAlarmRateAcrossBlocks: the average false alarm rate across
%               blocks of the experiment within each file
%           correctRejectionRateAcrossBlocks: the average hit rate across
%               blocks of the experiment within each file
%   
% Note that for early stereoscope alignment .log files the reported
% accuracy in the .log file may be inaccurate 9due to bug. This function
% recalculates accuracies from raw data and will produce the correct
% statistic.

% check inputs: 
if nargin<1
    error('AK_SSmechERP_logfileParse requires a subject''s initials as a string as an input')
end
if nargin<2
    data_dir = 'C:\Users\Alex Kale\Documents\MATLAB\CustomMurrayLabScripts\SSmechERP_behaviorAnalysis\data'; % default
end
% are they stings?
if ischar(subj_initials) == 0; 
    disp('subject intials are not a string')
end
if ischar(data_dir) == 0; 
    disp('data directory is not a string')
end

% % temporary inputs:
% subj_initials = 'vyc';
% data_dir = 'C:\Users\Alex Kale\Documents\MATLAB\CustomMurrayLabScripts\SSmechERP_behaviorAnalysis\data';

%% initialize variables for parsing logfiles using textscan

% parsing info
delimiter = '\t';
startRow = 4;
endRow = inf;
% format string for each line of text
formatSpec = '%*s%f%s%s%f%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% for classifying which experiment each logfile belongs to


%% find all files for this participant; categorize and load them

% find all logfiles for this subject
log_list = dir(fullfile(data_dir,[subj_initials '*.log']));
log_list = AK_sortStruct(log_list,1,'SSmechStereoAlign'); % sort files so that stereoAlign files are are in order
% preallocate output data structures and counters
stereoAlign = struct; % structs
C1 = struct;
SSmech = struct;
saCount = 0; % counters
c1Count = 0;
ssmCount = 0;
% cycle through logfiles
for iL = 1:length(log_list)

    % import events and times from tab-delimited log file:
    clear filename fileID dataArray dataArrayBlock data
    % generate filename
    filename = fullfile(data_dir,log_list(iL).name);
    % open the text file
    fileID = fopen(filename,'r');
    % read columns of data according to format string
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
    % close the text file.
    fclose(fileID);
    % organize data of interest for further processing
    dataArray([1, 4]) = cellfun(@(x) num2cell(x), dataArray([1, 4]), 'UniformOutput', false);
    data = [dataArray{1:end-1}];

    % classify this logfile as belonging to SSmech, C1 or StereoAlign experiment
    % parse data differently for each experiment, each gets its own data structure
    if ~isempty(regexp(log_list(iL).name,'.StereoAlign.','once')) % Stereoscope Alignment task
        clear stimIdx attemptBounds accuracy
        % increment counter for files of this type
        saCount = saCount + 1;
        % create index for when stimuli were presented
        stimIdx = find(strcmp(data(:,2),'Picture') & (strcmp(data(:,3),'left') | strcmp(data(:,3),'right') | strcmp(data(:,3),'correct')));
        % identify the number of sttempts at the alignment task in this .log file
        attemptBounds = [1; find(strcmp(data(:,2),'Picture') & ~cellfun(@isempty,cellfun(@(x) regexp(x,'Accuracy.'),data(:,3),'UniformOutput',false)))];
        % preallocate
        accuracy = zeros(length(attemptBounds) - 1,1); 
        for iA = 1:length(attemptBounds) - 1
            clear currentAttemptTrials correct
            % identify trials within current attempt
            currentAttemptTrials = stimIdx(ismember(stimIdx,attemptBounds(iA):attemptBounds(iA + 1)));
            % for each trial document where participant was correct and incorrect
            correct = zeros(length(currentAttemptTrials),1); % preallocate
            for iT = 1:length(currentAttemptTrials)
                if (strcmp(data(currentAttemptTrials(iT),3),'correct') && strcmp(data(currentAttemptTrials(iT) + 1,3),'1')) ||... % correct identification of proper alignment, or...
                        ((strcmp(data(currentAttemptTrials(iT),3),'left') || strcmp(data(currentAttemptTrials(iT),3),'right')) &&...
                        strcmp(data(currentAttemptTrials(iT) + 1,3),'2')) % correct identification of when stimulus is misaligned
                    correct(iT) = 1;
                end
            end
            accuracy(iA) = mean(correct);
        end
        % store data in structure
        stereoAlign(saCount).filename = log_list(iL).name;
        stereoAlign(saCount).date = log_list(iL).date;
        stereoAlign(saCount).accuracy = accuracy;
    elseif ~isempty(regexp(log_list(iL).name,'.C1test.','once')) % C1 experiment
        clear stimIdx blockBounds blockAccuracies noResponse correct
        % increment counter for files of this type
        c1Count = c1Count + 1;
        % create index for when stimuli were presented
        stimIdx = find(strcmp(data(:,2),'Picture') & ~cellfun(@isempty,cellfun(@(x) regexp(x,'Pic2.'),data(:,3),'UniformOutput',false)));
        % create index of block boundaries
        blockBounds = [find(strcmp(data(:,2),'Picture') & ~cellfun(@isempty,cellfun(@(x) regexp(x,'Block.'),data(:,3),'UniformOutput',false))); length(data(:,1))];
        % preallocate
        blockAccuracies = zeros(length(blockBounds) - 1,1);
        noResponse = cell(0);
        correct = cell(0);
        for iB = 1:length(blockBounds) - 1
            clear currentBlockTrials
            % identify trials within current block
            currentBlockTrials = stimIdx(ismember(stimIdx,blockBounds(iB):blockBounds(iB + 1)));
            % for each trial document where participant was correct and
            % incorrect and where they didn't respond
            % (there will be a number of responses equal to the number of
            % stimulus presentations minus one since this is a one-back
            % task)
            correct{iB} = zeros(length(currentBlockTrials),1); % preallocate
            noResponse{iB} = zeros(length(currentBlockTrials),1);
            for iT = 1:length(currentBlockTrials) 
                clear respIdx
                if iT == 1 || iT == length(currentBlockTrials) 
                    % first trial will never have response because one-back, and last trial responses were not recorded due to coding error
                    noResponse{iB}(iT) = 1;
                    correct{iB}(iT) = nan;
                else
                    % find the response for this trial
                    respIdx = find(strcmp(data(:,2),'Response') & ismember(1:length(data(:,1)),currentBlockTrials(iT):currentBlockTrials(iT + 1))',1); % find first response within current trial 
                    % has subject correctly identified whether current trial is
                    % the same vs different orientation from the previous trial
                    % (this is why we start counting trials at trial 2)
                    if isempty(respIdx) % is there a response this trial?
                        noResponse{iB}(iT) = 1;
                    elseif (((~isempty(regexp(data(currentBlockTrials(iT),3),'.or1.', 'once')) && ~isempty(regexp(data(currentBlockTrials(iT - 1),3),'.or1.', 'once'))) ||... % both orientation 1 or...
                            (~isempty(regexp(data(currentBlockTrials(iT),3),'.or1.', 'once')) && ~isempty(regexp(data(currentBlockTrials(iT - 1),3),'.or1.', 'once')))) &&... % both orientation two and...
                            strcmp(data(respIdx,3),'1')) ||... % response is '1' (same) or...
                            (((~isempty(regexp(data(currentBlockTrials(iT),3),'.or1.', 'once')) && ~isempty(regexp(data(currentBlockTrials(iT - 1),3),'.or2.', 'once'))) ||... % both different orientations
                            (~isempty(regexp(data(currentBlockTrials(iT),3),'.or2.', 'once')) && ~isempty(regexp(data(currentBlockTrials(iT - 1),3),'.or1.', 'once')))) &&... % and...
                            strcmp(data(respIdx,3),'2')) % response is '2' (different)
                        correct{iB}(iT) = 1;
                    end
                end
            end
            blockAccuracies(iB) = nanmean(correct{iB});
        end
        % store data in structure
        C1(c1Count).filename = log_list(iL).name;
        C1(c1Count).date = log_list(iL).date;
        C1(c1Count).correct = correct;
        C1(c1Count).noResponse = noResponse;
        C1(c1Count).accuracy = blockAccuracies;
        C1(c1Count).accuracyAcrossBlocks = mean(blockAccuracies); % each block should have same number of trials
    else % file belongs to SSmech main experiment
        % increment counter for files of this type
        ssmCount = ssmCount + 1;
        % create index for when stimuli were presented
        stimIdx = find(strcmp(data(:,2),'Picture') & ~cellfun(@isempty,cellfun(@(x) regexp(x,'Pic2.'),data(:,3),'UniformOutput',false)));
        % create index of block boundaries
        blockBounds = find(strcmp(data(:,2),'Picture') & ~cellfun(@isempty,cellfun(@(x) regexp(x,'Block.'),data(:,3),'UniformOutput',false))); 
        % check whether interblock interval displayed after last block; if so no need to include end of file idx
        quitIdx = find(strcmp(data(:,2),'Quit'));
        if quitIdx - blockBounds(end) > 2
            blockBounds = [blockBounds; quitIdx];
        end
        % preallocate
        hitRate = zeros(length(blockBounds) - 1,1);
        missRate = zeros(length(blockBounds) - 1,1);
        falseAlarmRate = zeros(length(blockBounds) - 1,1);
        correctRejectionRate = zeros(length(blockBounds) - 1,1);
        hits = cell(0);
        misses = cell(0);
        falseAlarms = cell(0);
        correctRejections = cell(0);
        isSignal = cell(0);
        for iB = 1:length(blockBounds) - 1
            clear currentBlockTrials
            % identify trials within current block
            currentBlockTrials = [stimIdx(ismember(stimIdx,blockBounds(iB):blockBounds(iB + 1))); blockBounds(iB + 1)]; % add end of block in order to catch last response
            % for each trial document hits, misses, false alarms, and
            % correct rejections
            hits{iB} = zeros(length(currentBlockTrials) - 1,1); % preallocate
            misses{iB} = zeros(length(currentBlockTrials) - 1,1);
            falseAlarms{iB} = zeros(length(currentBlockTrials) - 1,1);
            correctRejections{iB} = zeros(length(currentBlockTrials) - 1,1);
            isSignal{iB} = zeros(length(currentBlockTrials) - 1,1);
            for iT = 1:length(currentBlockTrials) - 1
                clear respIdx
                % find the response for this trial
                respIdx = find(strcmp(data(:,2),'Response') & ismember(1:length(data(:,1)),currentBlockTrials(iT):currentBlockTrials(iT + 1))',1); % find first response within current trial 
                % get the port code for this trial
                portCodeStrBounds = regexp(data{currentBlockTrials(iT),3},'[:]'); % the only ':'s in the event code will flank the port code 
                portCode = str2double(data{currentBlockTrials(iT),3}(portCodeStrBounds(1) + 1:portCodeStrBounds(2) - 1));
                % score plaid detection task
                if isempty(respIdx) % is there a response this trial? no:
                    if portCode > 99 % is the stimulus a plaid? yes:
                        % document plaid trial
                        isSignal{iB}(iT) = 1;
                        % document miss
                        misses{iB}(iT) = 1;
                    else % no:
                        % document correct rejection
                         correctRejections{iB}(iT) = 1;
                    end
                else % yes:
                    if portCode > 99 % is the stimulus a plaid? yes:
                        % document plaid trial
                        isSignal{iB}(iT) = 1;
                        % document hit
                         hits{iB}(iT) = 1;
                    else % no:
                        % document false alarm
                         falseAlarms{iB}(iT) = 1;
                    end
                end
            end
            hitRate(iB) = sum(hits{iB})/sum(isSignal{iB});
            missRate(iB) = sum(misses{iB})/sum(isSignal{iB});
            falseAlarmRate(iB) = sum(falseAlarms{iB})/sum(~isSignal{iB});
            correctRejectionRate(iB) = sum(correctRejections{iB})/sum(~isSignal{iB});
        end
        % store data in structure
        SSmech(ssmCount).filename = log_list(iL).name;
        SSmech(ssmCount).date = log_list(iL).date;
        SSmech(ssmCount).hits = hits;
        SSmech(ssmCount).misses = misses;
        SSmech(ssmCount).falseAlarms = falseAlarms;
        SSmech(ssmCount).correctRejections = correctRejections;
        SSmech(ssmCount).isSignal = isSignal;
        SSmech(ssmCount).hitRate = hitRate;
        SSmech(ssmCount).missRate = missRate;
        SSmech(ssmCount).falseAlarmRate = falseAlarmRate;
        SSmech(ssmCount).correctRejectionRate = correctRejectionRate;
        SSmech(ssmCount).hitRateAcrossBlocks = nanmean(hitRate);
        SSmech(ssmCount).missRateAcrossBlocks = nanmean(missRate);
        SSmech(ssmCount).falseAlarmRateAcrossBlocks = nanmean(falseAlarmRate);
        SSmech(ssmCount).correctRejectionRateAcrossBlocks = nanmean(correctRejectionRate);
    end
end

% aggregate statistics:
% the last measure accuracy attained per stereo alignment run 
for i = 1:length(stereoAlign)
    if ~isempty(stereoAlign(i).accuracy)
        overall.finalAccuracyStereoAlign(i) = vertcat(stereoAlign(i).accuracy(end));
    else
        overall.finalAccuracyStereoAlign(i) = nan;
    end
end
% accuracy for C1 task
overall.avgAccuracyC1 = mean(vertcat(C1.accuracyAcrossBlocks));
% measures for SSmech
overall.avgHitRateSSmech = nanmean(vertcat(SSmech.hitRateAcrossBlocks));
overall.avgMissRateSSmech = nanmean(vertcat(SSmech.missRateAcrossBlocks));
overall.avgFalseAlarmRateSSmech = nanmean(vertcat(SSmech.falseAlarmRateAcrossBlocks));
overall.avgCorrectRejectionRateSSmech = nanmean(vertcat(SSmech.correctRejectionRateAcrossBlocks));

end

