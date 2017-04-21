% ftap_ResponseVariability
% MurrayLab_2016
% Created by AMK 10/7/16

%% establish directory by subjects and fMRI sessions

top_dir = 'L:\MurrayLab\ASD\Data'; % set up base directory
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab\ftap_ResponseVariability'; % set directory for saving data and figures

% designate folders for invidual participants 
% subj_dirs = {'G101','G102','G103','G104','G105','G106','G107','G109','G110','G111','G112','G114','G117','G119','G307','G310','G311','G312','G313','G314','G315','G316','G317','G318','G319','G320','G322','G324','G326','G327','G328','G329','G332','G338'};
% get all subjects
subj_list = dir(fullfile(top_dir,'G*'));
subj_dirs = {subj_list.name};

conds = {'standard', 'variable'}; % create list of conditions for figure creation

%% load qualityTable.mat and filter out dqSubjects

load(fullfile(top_dir,'EyetrackingQuality.mat'));
subj_dirs = setdiff(subj_dirs,dqSubjects);

%% parse log files and store variability summary stats in structure

% preallocate
ftap = struct;

for iS = 1:length(subj_dirs); % cycle through subjects
    % generate directory in use
    use_dir = fullfile(top_dir,subj_dirs{iS},'ftap\log_files'); 
    if exist(use_dir,'dir')
        % generate structure containing information about log files
        log_list = dir(fullfile(use_dir,'*.log')); 
        log_list = AK_sortStruct(log_list,1,[subj_dirs{iS} '_']);
        if ~isempty(log_list)
            % preallocate
            stdRespIntervals = cell(1);
            varRespTimes = cell(1);
            for iL = 1:length(log_list) % cycle through log files
                % load ftap behavioral data
                clear timepts
                mat_name = [fullfile(use_dir,log_list(iL).name(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
                if ~exist(mat_name,'file') % check whether or not function output has already been saved
                    disp(['parsing ' fullfile(use_dir,log_list(iL).name)]) % message
                    [accuracy,resptimes,timepts] = AK_GABAinASD_ftap_logfileParse(fullfile(use_dir,log_list(iL).name)); % run log file through parsing function
                    save(mat_name,'accuracy','resptimes','timepts'); % save function output
                else
                    disp(['loading ' mat_name]) % message
                    load(mat_name,'timepts') % load file containing timepts array
                end 
                % calculate timing variability for 'standard' (predictable) trials
                clear stdRespIdx stdRespTimestamps
                stdRespIdx = strcmp(timepts(:,5),'hitStd'); % index for rows of timepts that correspond to responses to the standard cond
                stdRespTimestamps = vertcat(timepts{stdRespIdx,3}); % raw timestamps
                stdRespIntervals{iL} = (stdRespTimestamps(2:end)-stdRespTimestamps(1:end-1)); % differences between consecutive timestamps
                stdRespIntervals{iL}(stdRespIntervals{iL}>30325) = []; % where gap between responses is more than a second larger than the gap between stim presentations, intervals correspond to block gaps and should be removed from analysis
                stdRespIntervals{iL} = stdRespIntervals{iL}./10000; % convert to seconds
                % calculate response time variability for 'variable' (random) trials
                clear varRespIdx
                varRespIdx = strcmp(timepts(:,5),'hitVar'); % index for rows of timepts that correspond to responses to the variable cond
                varRespTimes{iL} = vertcat(timepts{varRespIdx,4}); % response times calculated in parsing function 
                varRespTimes{iL} = varRespTimes{iL}./10000; % convert to seconds
            end
            % concatenate across logfiles
            stdResp = vertcat(stdRespIntervals{:});
            varResp = vertcat(varRespTimes{:});
            % store summary stats in structure
            ftap(iS).subj = subj_dirs{iS};
            ftap(iS).predictableM = nanmean(stdResp);
            ftap(iS).predictableSD = nanstd(stdResp);
            ftap(iS).randomM = nanmean(varResp);
            ftap(iS).randomSD = nanstd(varResp);
        else
            % put subject code in structure
            ftap(iS).subj = subj_dirs{iS};
        end
    else
        % put subject code in structure
        ftap(iS).subj = subj_dirs{iS};
    end
end

%% prepare desired output format

% column key
columnKey = {'predictableM','predictableSD','randomM','randomSD'};

for i = 1:length(ftap)
    % row key
    rowKey{i} = ftap(i).subj;
    % means and standard deviations
    if isempty(ftap(i).predictableM)
        ftapVariability(i,1:2) = nan;
    else
        ftapVariability(i,1) = ftap(i).predictableM;
        ftapVariability(i,2) = ftap(i).predictableSD;
    end
    if isempty(ftap(i).randomM)
        ftapVariability(i,3:4) = nan;
    else
        ftapVariability(i,3) = ftap(i).randomM;
        ftapVariability(i,4) = ftap(i).randomSD;
    end
end

%% save

save(fullfile(home_dir,'ftapVariability.mat'),'ftapVariability','rowKey','columnKey');
