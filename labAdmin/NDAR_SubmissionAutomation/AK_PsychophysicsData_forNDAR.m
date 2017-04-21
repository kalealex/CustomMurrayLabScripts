function AK_PsychophysicsData_forNDAR( data_filename, subjInfo_filename )
%AK_PsychophysicsData_forNDAR prepares a .csv file containing appropriately 
% formatted data for all psychophysics sessions for which data upload is 
% necessary and saves this file to the NDAR submission folder with the 
% fullfile name specified in the first function input (data_filename). It 
% does this for the subjects whose information is contained in the file 
% designated in the second function input (subjInfo_filename). See 
% AK_prepareNDARsubmission.m documentation for more detailed instructions 
% on the NDAR submission process and file formatting.
%
% HISTORY:
% used to be Psychophysics_Data_for_NDAR
% Murray Lab 2015
% Created by Alex Kale 1/8/16
% Edited 6/15/16

%% establish directory and grab subject info

top_dir = 'L:\MurrayLab\ASD\Data'; % set up base directory
psychophys_dirs = {'Psychophysics','Psychophysics2'}; % subfolders for psychophysics

% subjInfo_filename = fullfile(top_dir,'NDARnotes','NDAR_SubjectInfo.xlsx'); % create full filename for NDAR subject info spreadsheet
[~,~,subjInfo] = xlsread(subjInfo_filename); % load subject info
subj_dirs = subjInfo(2:end,1); % designate subject ids

% nMSconds = vertcat(subjInfo{2:end,5}); % load number of conditions in motion staircase
dateFormat = 'mm/dd/yyyy'; % string for date formatting

%% create cell array for data output

NDARdata = cell(6*length(subj_dirs)*196+1,21); % rows: 6 .mat files * #subjects * 196 possible conditions per file
NDARdata(1,1:2) = {'psef','1'};
NDARdata(2,:) = {'subjectkey','src_subject_id','interview_age','interview_date','gender','experiment_name','experiment_notes','rundata_roomno','rundata_responsetime','rundata_conditionorder','rundata_sideorder','savestair_stimrange','savestair_guess','savestair_lapse','savestair_x','savestair_response','savestair_threshold','savestair_slope','savestair_sethreshold','savestair_seslope','response_number'};

%% load data and organize in table, per file, per subject

NDARrow = 3; % designate variable for indexing NDARdata rows
nFiles = zeros(length(subj_dirs),1); % create vector for storing number of mat files per subject
nSavestair = cell(length(subj_dirs),1); % create cell array to document number of elements in stairsave, per file, per subject

for iS = 1:length(subj_dirs); % cycle through subjects
    subj_dir = fullfile(top_dir,subj_dirs{iS});
    if exist(subj_dir,'dir')==0
        disp(['could not find subject folder for' subj_dirs{iS}]) % message
    elseif exist(subj_dir,'dir')==7
        for iPSY = 1:length(psychophys_dirs); % cycle through potential psychophysics folders
            use_dir = fullfile(subj_dir,psychophys_dirs{iPSY}); % generate directory in use
            if exist(use_dir,'dir')==0 % psychophysics directory doesn't exist
                disp([use_dir ' is not an existing directory']); % message
            else 
                mat_list = dir(fullfile(use_dir,'*.mat')); % generate structure containing information about mat files
                if ~isempty(mat_list)
                    clear MotionStairStr ContrastStairStr rawdata_mat_name runData savestair
                    MotionStairStr = [subj_dirs{iS} '_motion_staircase*']; % create string for .mat file of interest
                    ContrastStairStr = [subj_dirs{iS} '_contrast_staircase*']; % create string for .mat file of interest
                    for iM = 1:length(mat_list); % cycle through raw data files
                        if ~isempty(regexp(mat_list(iM).name,regexptranslate('wildcard',MotionStairStr), 'once'))
                            rawdata_mat_name = fullfile(use_dir,mat_list(iM).name); % create string name of raw data .mat file
                            disp(['loading ' rawdata_mat_name]) % message
                            load(rawdata_mat_name,'runData','savestair') % load file containing raw psychophysics data
                            % create index for conditions to extract responses from
                            % savestair
                            condIndex = cell(length(savestair),1);
                            for iMSC = 1:length(savestair); % cycle through motion staircase conditions
                                condIndex{iMSC} = find(runData.conditionOrder==iMSC);
                            end
%                             condIndex = cell(nMSconds(iS),1);
%                             for iMSC = 1:nMSconds(iS); % cycle through motion staircase conditions
%                                 condIndex{iMSC} = find(runData.conditionOrder==iMSC);
%                             end

                            for C = 1:length(condIndex); % cycle through conditions
                                for NP = 1:length(condIndex{C}); % cycle through number of presentations of each condition
                                    if C == length(savestair) % for catch trials only
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,16} = savestair(C).response(NP); % store response
                                    else
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,12} = savestair(C).stimRange(NP); % store stimRange
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,13} = savestair(C).guess(NP); % store guess
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,14} = savestair(C).lapse(NP); % store lapse
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,15} = savestair(C).x(NP); % store x
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,16} = savestair(C).response(NP); % store response
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,17} = savestair(C).threshold(NP); % store threshold
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,18} = savestair(C).slope(NP); % store slope
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,19} = savestair(C).seThreshold(NP); % store seThreshold
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,20} = savestair(C).seSlope(NP); % store seSlope
                                    end
                                end
                            end

                            for iT = 1:length(runData.conditionOrder) % cycle through trials
                                NDARdata{NDARrow,1} = subjInfo{iS+1,2}; % store GUID
                                NDARdata{NDARrow,2} = subj_dirs{iS}; % store subject id
                                NDARdata{NDARrow,3} = AK_calculateInterviewAge_forNDAR(subjInfo{iS+1,3}, datestr(mat_list(iM).date,dateFormat)); % store age
                                NDARdata{NDARrow,5} = subjInfo{iS+1,4}; % store gender

                                NDARdata{NDARrow,4} = datestr(mat_list(iM).date,dateFormat); % store interview date
                                NDARdata{NDARrow,6} = mat_list(iM).name(6:end-4); % store experiment name
                                NDARdata{NDARrow,8} = runData.roomNo; % store room number
                                NDARdata{NDARrow,9} = runData.responseTime(iT); % store response times
                                NDARdata{NDARrow,10} = runData.conditionOrder(iT); % store condition order
                                NDARdata{NDARrow,11} = runData.sideOrder(iT); % store side order

                                NDARdata{NDARrow,21} = iT; % store response number

                                NDARrow = NDARrow+1; % move to next row of NDARdata
                            end

                            nFiles(iS) = nFiles(iS)+1; % document number of files per subject
                            nSavestair{iS}(iM) = length(savestair); % document number of staircases
                        elseif ~isempty(regexp(mat_list(iM).name,regexptranslate('wildcard',ContrastStairStr), 'once'))
                            rawdata_mat_name = fullfile(use_dir,mat_list(iM).name); % create string name of raw data .mat file
                            disp(['loading ' rawdata_mat_name]) % message
                            load(rawdata_mat_name,'runData','savestair') % load file containing raw psychophysics data
                            % create index for conditions to extract responses from
                            % savestair
                            condIndex = cell(length(savestair),1);
                            for iCDC = 1:length(savestair); % cycle through condrast detection conditions
                                condIndex{iCDC} = find(runData.conditionOrder==iCDC);
                            end
%                             condIndex = cell(3,1);
%                             condIndex{1} = find(runData.conditionOrder==1);
%                             condIndex{2} = find(runData.conditionOrder==2);
%                             condIndex{3} = find(runData.conditionOrder==3);

                            for C = 1:length(condIndex); % cycle through conditions
                                for NP = 1:length(condIndex{C}); % cycle through number of presentations of each condition
                                    if C == length(savestair) % for catch trials only
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,16} = savestair(C).response(NP); % store response
                                    else
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,12} = savestair(C).stimRange(NP); % store stimRange
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,13} = savestair(C).guess(NP); % store guess
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,14} = savestair(C).lapse(NP); % store lapse
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,15} = savestair(C).x(NP); % store x
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,16} = savestair(C).response(NP); % store response
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,17} = savestair(C).threshold(NP); % store threshold
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,18} = savestair(C).slope(NP); % store slope
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,19} = savestair(C).seThreshold(NP); % store seThreshold
                                        NDARdata{NDARrow+condIndex{C}(NP)-1,20} = savestair(C).seSlope(NP); % store seSlope
                                    end
                                end
                            end

                            for iT = 1:length(runData.conditionOrder) % cycle through trials
                                NDARdata{NDARrow,1} = subjInfo{iS+1,2}; % store GUID
                                NDARdata{NDARrow,2} = subj_dirs{iS}; % store subject id
                                NDARdata{NDARrow,3} = AK_calculateInterviewAge_forNDAR(subjInfo{iS+1,3}, datestr(mat_list(iM).date,dateFormat)); % store age
                                NDARdata{NDARrow,5} = subjInfo{iS+1,4}; % store gender

                                NDARdata{NDARrow,4} = datestr(mat_list(iM).date,dateFormat); % store interview date
                                NDARdata{NDARrow,6} = mat_list(iM).name(6:end-4); % store experiment name
                                NDARdata{NDARrow,8} = runData.roomNo; % store room number
                                NDARdata{NDARrow,9} = runData.responseTime(iT); % store response times
                                NDARdata{NDARrow,10} = runData.conditionOrder(iT); % store condition order
                                NDARdata{NDARrow,11} = runData.orientOrder(iT); % store side order (orient order)

                                NDARdata{NDARrow,21} = iT; % store response number

                                NDARrow = NDARrow+1; % move to next row of NDARdata
                            end

                            nFiles(iS) = nFiles(iS)+1; % document number of files per subject
                            nSavestair{iS}(iM) = length(savestair); % document number of staircases
                        else
                            disp(['failure to recognize file' mat_list(iM).name]) % message
                        end
                    end
                else
                    disp(['mat_list is empty for ' psychophys_dirs{iPSY}]) % message
                end
            end
        end
    end
end

%% find and delete blank rows

isEmptyIndex = zeros(length(NDARdata(:,1)),1);

for iNDAR = 1:length(NDARdata(:,1)); % cycle through rows of NDARdata
    if mean(cellfun(@isempty,NDARdata(iNDAR,:))) == 1; % is the entire row empty?
        isEmptyIndex(iNDAR) = 1;
    end
end

isEmptyIndex = logical(isEmptyIndex); % make logical

NDARdata(isEmptyIndex,:) = []; % clear empty rows

%% write to csv file

cell2csv(data_filename,NDARdata);


end

