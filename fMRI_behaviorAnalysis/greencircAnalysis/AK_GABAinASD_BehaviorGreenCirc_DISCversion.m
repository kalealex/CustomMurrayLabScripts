function [ data ] = AK_GABAinASD_BehaviorGreenCirc_DISCversion( subjectCode, scanNumber, experimentDirectory )
%AK_GABAinASD_BehaviorGreenCirc arguments are subject code as a string, 
%scanNumber, and an optional argument designating the directory where
%the experiment files are located (the function should default to the
%right one). The function outputs behavioral data for the green circle task
%for the 'GABA in ASD' project. 
%   Input: 
%    subject code as a string 
%    scan number as a double (same one entered into presentnation)
%    experiment directory [optional], otherwise function will default to 
%    ssMotion folder
%   The output 'data' structure is organized with one set of fields
%   per subject(i.e. data(n).subject = subjectCodeCellStr{n}). Within each
%   field (except for 'subject' and 'fMRIsession') there are #fMRIsessions 
%   by #logfiles dimensions (i.e. data(n).accuracy(2,1) = table of accuracy
%   info for first logfile of second fMRI session). For each logfile, there
%   are fields for logfile name, condition names, a pulse count, a table of 
%   accuracy stats separated by condition, and a list of response times
%   separated by condition. Accuracy and response time data is also 
%   provided in the aggregate under the condition named 'all'.


% Check arguments
if nargin < 2;
    error('AK_GABAinASD_BehaviorGreenCirc_DISCversion requires subject code (string) and scan number (double) as arguments');
end
if nargin < 3;
   experimentDirectory = [];
end

% establish directory and logfile finding string 
if ~exist('experimentDirectory','var') || isempty(experimentDirectory);
    experimentDirectory = 'enter experiment directory here';
end
use_dir = fullfile(experimentDirectory,'logfiles'); % set up logfile directory
scanStr = [subjectCode '_' num2str(scanNumber)]; % create string to be used to find logfile

% preallocate for speed
data = struct;

% find and parse log file; store in data structure
clear log_list
log_list = dir(fullfile(use_dir,[scanStr '*.log'])); % generate structure containing information about log files
if ~isempty(log_list)
    for iL = 1:length(log_list) % cycle through log files
        % parse logfile
        clear accuracy resptimes timepts pulseCount
        disp(['parsing ' fullfile(use_dir,log_list(iL).name)]) % message
        [accuracy,resptimes,timepts,pulseCount] = AK_GABAinASD_logfileParse(fullfile(use_dir,log_list(iL).name)); % run log file through parsing function

        % directory information    
        data.subject = subjectCode;       
        data.logFile{iL} = log_list(iL).name; % store name of logfile

        % create indices for determining accuracy information by
        % condition
        clear hit_index miss_index false_alarm_index greencirc_index other_index resp_index cond_list
        hit_index = find(strcmp(timepts(2:end,5),'hit')); % create index of hits
        miss_index = find(strcmp(timepts(2:end,5),'miss')); % create index of misses
        false_alarm_index = find(strcmp(timepts(2:end,5),'false alarm')); % create index of false alarms
        greencirc_index = find(strcmp(timepts(2:end,2),'greencirc')); % create index of green circle presentations
        other_index = find(strcmp(timepts(2:end,2),'other')); % create index of other shapes
        resp_index = find(strcmp(timepts(2:end,2),'Response')); % create index of responses
        nonempty_conds = timepts(~cellfun(@isempty,timepts(:,6)),6); % create list of all nonempty conditions including label
        cond_list = unique(nonempty_conds(2:end)); % create a list of unique condition labels; exclude label in row 1

        % preallocate for speed
        clear cond_index cond_respTimes cond_hits cond_misses cond_false_alarms cond_greenCircs cond_otherShapes cond_hit_rate cond_miss_rate cond_false_alarm_rate cond_dPrime 
        cond_index = cell(length(cond_list),1);
        cond_respTimes = cell(length(cond_list),1);
        cond_hits = zeros(length(cond_list),1);
        cond_misses = zeros(length(cond_list),1);
        cond_false_alarms = zeros(length(cond_list),1);
        cond_greenCircs = zeros(length(cond_list)+1,1);
        cond_otherShapes = zeros(length(cond_list)+1,1);
        cond_hit_rate = zeros(length(cond_list)+1,1);
        cond_miss_rate = zeros(length(cond_list)+1,1);
        cond_false_alarm_rate = zeros(length(cond_list)+1,1);
        cond_dPrime = zeros(length(cond_list)+1,1);
        for iC = 1:length(cond_list); % cycle through unique conditions
            cond_index{iC} = find(strcmp(timepts(2:end,6),cond_list(iC))); % create index cells of timepts(2:end,:) within each condition
            % separate response times by condition
            if ~isempty(resp_index) && ~isempty(intersect(resp_index,cond_index{iC})); % check for responses in cond iC
                cond_respTimes{iC} = timepts{intersect(resp_index,cond_index{iC})+1,4};
            else
                cond_respTimes{iC} = [];
            end
            % count hits, misses, false alarms, and # stimulus
            % presentations (signal and not) per condition
            cond_hits(iC) = length(intersect(hit_index,cond_index{iC}));
            cond_misses(iC) = length(intersect(miss_index,cond_index{iC}));
            cond_false_alarms(iC) = length(intersect(false_alarm_index,cond_index{iC}));
            cond_greenCircs(iC) = length(intersect(greencirc_index,cond_index{iC}));
            cond_otherShapes(iC) = length(intersect(other_index,cond_index{iC}));
            % determine hit rate, miss rate, false alarm rate, and d prime per condition
            cond_hit_rate(iC) = cond_hits(iC)/cond_greenCircs(iC);
            cond_miss_rate(iC) = cond_misses(iC)/cond_greenCircs(iC);
            cond_false_alarm_rate(iC) = cond_false_alarms(iC)/cond_otherShapes(iC);
            cond_dPrime(iC) = AK_dPrime(cond_greenCircs(iC),cond_otherShapes(iC),cond_hits(iC),cond_false_alarms(iC));
        end

        % add data for all trials across conditions
        cond_greenCircs(iC+1) = length(greencirc_index);
        cond_otherShapes(iC+1) = length(other_index);
        cond_hit_rate(iC+1) = length(hit_index)/length(greencirc_index);
        cond_miss_rate(iC+1) = length(miss_index)/length(greencirc_index);
        cond_false_alarm_rate(iC+1) = length(false_alarm_index)/length(other_index);
        cond_dPrime(iC+1) = AK_dPrime(length(greencirc_index),length(other_index),length(hit_index),length(false_alarm_index));

        % designate labels for accuracy table
        clear accuracyTable_condLabels accuracyTable_statLabels
        accuracyTable_condLabels = [cond_list;{'all'}]';
        accuracyTable_statLabels = {'condition','hitRate','missRate','falseAlarmRate','nGreenCircles','nOtherShapes'};

        % document table of accuracy statistics
        data.accuracy{iL} = table(cond_dPrime,cond_hit_rate,cond_miss_rate,cond_false_alarm_rate,cond_greenCircs,cond_otherShapes,'RowNames',accuracyTable_condLabels,'VariableNames',accuracyTable_statLabels);
        data.accuracy{iL} = cell(length(accuracyTable_condLabels)+1,length(accuracyTable_statLabels));
        data.accuracy{iL}(1,:) = accuracyTable_statLabels;
        data.accuracy{iL}(2:end,1) = accuracyTable_condLabels;
        data.accuracy{iL}(2:end,2) = arrayfun(@num2cell,cond_hit_rate);
        data.accuracy{iL}(2:end,3) = arrayfun(@num2cell,cond_miss_rate);
        data.accuracy{iL}(2:end,4) = arrayfun(@num2cell,cond_false_alarm_rate);
        data.accuracy{iL}(2:end,5) = arrayfun(@num2cell,cond_greenCircs);
        data.accuracy{iL}(2:end,6) = arrayfun(@num2cell,cond_otherShapes);

        % document response times per condition in structure
        data.conds{iL} = [cond_list;{'all'}]';
        data.respTimes{iL} = [cond_respTimes;{resptimes.all}]';

        % document pulse count for each logfile as a means of
        % detecting incomplete logfiles
        data.pulseCount(iL) = pulseCount; 
    end
else
    disp(['Logfile does not exist in ' use_dir]); % message
end


end
