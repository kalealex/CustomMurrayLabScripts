function [ data ] = AK_GABAinASD_BehaviorGreenCirc( subjectCodeCellStr, logfileDirectory, homeDirectory )
%AK_GABAinASD_BehaviorGreenCirc arguments are subject code[s] as a cell array of
%strings, as well as optional arguments designating the directory where
%logfiles are located and the directory where the data structure for the
%designated subjects should be saved. The function outputs behavioral data
%for the green circle task for the 'GABA in ASD' project. 
%   Input subject code[s] as cell array of strings. 
%   Input logfile directory [optional], otherwise function will default to L
%   drive. 
%   Input home directory [optional] if you want the function to save output
%   automatically. 
%   The output 'data' structure is organized with one set of fields
%   per subject(i.e. data(n).subject = subjectCodeCellStr{n}). Within each
%   field (except for 'subject' and 'fMRIsession') there are #fMRIsessions 
%   by #logfiles dimensions (i.e. data(n).accuracy(2,1) = table of accuracy
%   info for first logfile of second fMRI session). For each logfile, there
%   are fields for logfile name, condition names, a pulse count, a table of accuracy stats
%   separated by condition, and a list of response times separated by
%   condition. Accuracy and response time data is also provided in the
%   aggregate under the condition named 'all'.

% Check arguments
if nargin < 1;
    error('AK_GABAinASD_BehaviorGreenCirc requires subject code[s] (cell array of strings) as an argument');
end
if nargin < 3;
   homeDirectory = [];
end
if nargin < 2;
   logfileDirectory = [];
end


% establish directory by subjects and fMRI sessions
if ~exist('logfileDirectory','var') || isempty(logfileDirectory);
    logfileDirectory = 'L:\MurrayLab\ASD\Data';
end
top_dir = logfileDirectory; % set up base directory
if exist('homeDirectory','var') || ~isempty(homeDirectory);
    home_dir =  homeDirectory; % set directory for saving data
end
subj_dirs = subjectCodeCellStr; % designate folders for invidual participants 
fmri_dirs = {'fMRI1','fMRI2'}; % designate subfolders

% preallocate for speed
data(1:length(subj_dirs)) = struct;

% parse log files and store in data structure
for iS = 1:length(subj_dirs); % cycle through subjects
    for iF = 1:length(fmri_dirs); % cycle through fmri sessions
        clear use_dir log_list
        use_dir = fullfile(top_dir,subj_dirs{iS},fmri_dirs{iF},'log_files'); % generate directory in use
        log_list = dir(fullfile(use_dir,'*.log')); % generate structure containing information about log files
        log_list = AK_sortStruct(log_list,1,[subj_dirs{iS} '_']); % sort structure by scan# in ascending order
        if ~isempty(log_list)
            for iL = 1:length(log_list) % cycle through log files
                clear accuracy resptimes timepts pulseCount
                mat_name = [fullfile(use_dir,log_list(iL).name(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
                    if ~exist(mat_name,'file') % check whether or not function output has already been saved
                        disp(['parsing ' fullfile(use_dir,log_list(iL).name)]) % message
                        [accuracy,resptimes,timepts,pulseCount] = AK_GABAinASD_logfileParse(fullfile(use_dir,log_list(iL).name)); % run log file through parsing function
                        save(mat_name,'accuracy','resptimes','timepts','pulseCount'); % save function output
                    else
                        disp(['loading ' mat_name]) % message
                        load(mat_name,'accuracy','resptimes','timepts','pulseCount') % load file containing timepts array
                    end
                
                % directory information    
                data(iS).subject = subj_dirs{iS};    
                data(iS).fMRIsession{iF,1} = fmri_dirs{iF};    
                data(iS).logFile{iF,iL} = log_list(iL).name; % store name of logfile
                
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
                accuracyTable_statLabels = {'dPrime','hitRate','missRate','falseAlarmRate','nGreenCircles','nOtherShapes'};
                
                % document table of accuracy statistics
                data(iS).accuracy{iF,iL} = table(cond_dPrime,cond_hit_rate,cond_miss_rate,cond_false_alarm_rate,cond_greenCircs,cond_otherShapes,'RowNames',accuracyTable_condLabels,'VariableNames',accuracyTable_statLabels);
                
                % document response times per condition in structure
                data(iS).conds{iF,iL} = [cond_list;{'all'}]';
                data(iS).respTimes{iF,iL} = [cond_respTimes;{resptimes.all}]';
                
                % document pulse count for each logfile as a means of
                % detecting incomplete logfiles
                data(iS).pulseCount(iF,iL) = pulseCount; 
            end
        else
            disp(['Logfiles do not exist in ' use_dir]); % message
        end
    end
end

% save data structure
if exist(homeDirectory,'var');
    if exist(homeDirectory,'file');
        data_dir = home_dir; % set directory for saving data
        data_name = fullfile(data_dir,'fMRI_Behavioral_Data.mat'); % create string to name .mat file which saves data
        if ~exist(data_name,'file') % check whether or data have already been saved
            disp(['saving ' data_name]) % message
            save(data_name,'data'); % save data
        end
    end
end

end
