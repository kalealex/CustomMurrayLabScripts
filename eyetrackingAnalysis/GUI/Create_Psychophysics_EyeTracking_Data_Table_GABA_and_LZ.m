% Create_Psychophysics_EyeTracking_Data_Table_GABA_and_LZ
% MurrayLab_2017
% Created by AMK 3/20/17

% switches:
%create table for good quality data only?
onlyGood = true;
%only analyze the data (false to create table)
analyzeOnly = false;

%% designate directories

top_dir = 'L:\MurrayLab';
gaba_dir = 'ASD\Data';
lz_dir = 'Lorazepam\Data';
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab';
save_dir = 'L:\MurrayLab\DataTablesForGUI';

%% add paths fir AMK functions

analysis_path = fullfile(top_dir,gaba_dir,'Analysis_scripts','AMK_Functions');
addpath(genpath(analysis_path));

%% find subjects to include in analysis

% for GABA in ASD
gaba_subjects = dir(fullfile(top_dir,gaba_dir,'G*')); % structure of all matching items in directory
gaba_subjects = {gaba_subjects(:).name}; % cell array of subject codes/ folder names
gaba_subjects = [gaba_subjects {'KG122'}]; % add kiddo
gaba_subjects(cellfun(@(x) exist(x,'dir')==0,cellfun(@(x) fullfile(top_dir,gaba_dir,x),gaba_subjects,'UniformOutput',false))) = []; % remove elements which are not valid directories
% for Lorazepam
lz_subjects = dir(fullfile(top_dir,lz_dir,'*L*')); % structure of all matching items in directory
lz_subjects = {lz_subjects(:).name}; % cell array of subject codes/ folder names
lz_subjects(cellfun(@(x) exist(x,'dir')==0,cellfun(@(x) fullfile(top_dir,lz_dir,x),lz_subjects,'UniformOutput',false))) = []; % remove elements which are not valid directories

if ~analyzeOnly
    %% designate table structure and field names

    
    % name columns of eye tracking data tables (condition name: statistic name)
    gaba_columns = {'subject','run#','directory','behavioral data filename','eyetracking data filename'};
    lz_columns = [gaba_columns(1) {'session'} gaba_columns(2:end)];
    statList = {'no tracking time','fix time','drift corrected fix time','saccade and blink time','distance M','distance SD',...
        'angle M','angle SD','drift corrected distance M','drift corrected distance SD','drift corrected angle M','drift corrected angle SD',...
        'n saccades','n big saccades','n blinks','saccade amplitude M','saccade amplitude SD','saccade velocity M','saccade velocity SD'};
    motion_CondList = {'3% contrast, 2 degrees','3% contrast, 12 degrees','98% contrast, 2 degrees','98% contrast, 12 degrees','3% contrast, 1 degree','98% contrast, 1 degree','catch trials','all'};
    motion_CondListShort = {'3% contrast, 2 degrees','3% contrast, 12 degrees','98% contrast, 2 degrees','98% contrast, 12 degrees','catch trials','all'};
    contrast_gabaCondList = {'vertical','horizontal','catch trials','all'};
    contrast_lzCondList = {'staircase1','staircase2','staircase3','staircase4','catch trials','all'};

    % preallocate
    motion_columns = cell(1,length(statList)*length(motion_CondList));
    contrastGABA_columns = cell(1,length(statList)*length(contrast_gabaCondList));
    contrastLZ_columns = cell(1,length(statList)*length(contrast_lzCondList));
    % put together stats and conds for column headers i.e., <cond>: <stat>
    for iS = 1:length(statList)
        for iM = 1:length(motion_CondList)
            motion_columns{length(statList) * (iM - 1) + iS} = [motion_CondList{iM} ': ' statList{iS}];
        end
        for iCG = 1:length(contrast_gabaCondList)
            contrastGABA_columns{length(statList) * (iCG - 1) + iS} = [contrast_gabaCondList{iCG} ': ' statList{iS}];
        end
        for iCL = 1:length(contrast_lzCondList)
            contrastLZ_columns{length(statList) * (iCL - 1) + iS} = [contrast_lzCondList{iCL} ': ' statList{iS}];
        end
    end

    % create cell arrays to store data and save into spreadsheets
    % give each sheet 50 rows of buffer
    GABAmotion = cell(4 * length(gaba_subjects) + 50, length(gaba_columns) + length(motion_columns));
    LZmotion = cell(4 * length(lz_subjects) + 50, length(lz_columns) + length(motion_columns));
    GABAcontrast = cell(4 * length(gaba_subjects) + 50, length(gaba_columns) + length(contrastGABA_columns));
    LZcontrast = cell(4 * length(lz_subjects) + 50, length(lz_columns) + length(contrastLZ_columns));

    % add column names to cell arrays
    GABAmotion(1,:) = [gaba_columns,motion_columns];
    LZmotion(1,:) = [lz_columns,motion_columns];
    GABAcontrast(1,:) = [gaba_columns,contrastGABA_columns];
    LZcontrast(1,:) = [lz_columns,contrastLZ_columns];

    % create counters for rows of each cell array
    GABAmotionRow = 2;
    LZmotionRow = 2;
    GABAcontrastRow = 2;
    LZcontrastRow = 2;

    %% unblind Lorazepam sessions

    % load unblinding key
    cd(fullfile(top_dir,'Lorazepam/analysis_code'));
    load('blind_order.mat');
    % use blind_order key to index into a treatment list creating a cell array
    % of labeled Lorazepam sessions (drug vs placebo) by subject
    treatmentList = {'drug','placebo'};
    treatment = treatmentList(blind_order)';
end

%% pull data from data structure into appropriate spreadsheets

% build tables for GABA data
for iS = 1:length(gaba_subjects)
    % parse eye tracking data and save data structure for each subject individually
    clear eyetracking_data_file eyetracking_data_dir data
    eyetracking_data_file = [gaba_subjects{iS} '_Psychophysics_EyeTracking_Data.mat'];
    eyetracking_data_dir = fullfile(top_dir,gaba_dir,gaba_subjects{iS});
    cd(eyetracking_data_dir);
    if ~exist(eyetracking_data_file,'file') % look for saved data
        data = AK_MotionContrastPsychophysics_EvaluateFixation(gaba_subjects(iS));
        save(eyetracking_data_file,'data');
    elseif ~analyzeOnly
        disp(['loading ' eyetracking_data_file]) % message
        load(eyetracking_data_file,'data');
    end
    
    if ~analyzeOnly   
        % get motion data 
        expIdx = 1; % motion
        if ~isempty(data.experiment(expIdx).ascFilename) && (iscell(data.experiment(expIdx).ascFilename) || ~isnan(data.experiment(expIdx).ascFilename)) % is there eye tracking data for this experiment?
            for iRun = 1:length(data.experiment(expIdx).runName)
                if onlyGood
                    if data.experiment(expIdx).quality{iRun} == 1
                        % checks to advance row counter and avoid overwriting
                        if all(~cellfun(@isempty,GABAmotion(GABAmotionRow,length(gaba_columns)+1:length(GABAmotion(1,:))))) % check to advance counter
                            GABAmotionRow = GABAmotionRow+1; % advance counter
                        elseif any(~cellfun(@isempty,GABAmotion(GABAmotionRow,length(gaba_columns)+1:length(GABAmotion(1,:))))) % check that data is not being overwritten
                            disp(['Avoided overwriting for asc file: ' data.experiment(expIdx).ascFilename{iRun} ' at row ' num2str(GABAmotionRow) ' of GABAmotion']) % message
                            GABAmotionRow = GABAmotionRow+1; % advance counter
                        end

                        % store metadata about eye tracking data
                        GABAmotion{GABAmotionRow,1} = gaba_subjects{iS}; % store subject name
                        GABAmotion{GABAmotionRow,2} = iRun; % store run number
                        GABAmotion{GABAmotionRow,3} = data.directory; % store directory
                        GABAmotion{GABAmotionRow,4} = data.experiment(expIdx).behavioralDataFilename{iRun}; % store behavioral data filename
                        GABAmotion{GABAmotionRow,5} = data.experiment(expIdx).ascFilename{iRun}; % store eyetracking data filename

                        % store stats by condition
                        for iCond = 1:length(motion_CondList)
                            if (size(data.experiment(expIdx).runStats{iRun},1) == 6 && (iCond == 5 || iCond == 6)) || ~istable(data.experiment(expIdx).runStats{iRun})
                                % conditions 5 and 6 (1 degree stimuli were
                                % not run for early GABA participants) 
                            else
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 6) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'noTrackTime')); % store no tracking time for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 7) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'fixationTime')); % store fixation time for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 8) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedFixationTime')); % store drift-corrected fixation time for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 9) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'saccadeTime')); % store saccade and blink time for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 10) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'distanceMean')); % store distance M for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 11) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'distanceSD')); % store distance SD for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 12) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'angleMean')); % store angle M for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 13) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'angleSD')); % store angle SD for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 14) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedDistanceMean')); % store drift-corrected distance M for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 15) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedDistanceSD')); % store drift-corrected distance SD for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 16) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedAngleMean')); % store drift-corrected angle M for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 17) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedAngleSD')); % store drift-corrected angle SD for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 18) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeN')); % store number of saccades for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 19) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'bigSaccadeN')); % store number of big saccades for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 20) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'BlinkN')); % store number of saccades for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 21) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeAmplitudeMean')); % store saccade amplitude mean for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 22) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeAmplitudeSD')); % store saccade amplitude SD for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 23) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeVelocityMean')); % store saccade velocity mean for condition
                                GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 24) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeVelocitySD')); % store saccade velocity SD for condition
                            end
                        end    
                    end
                else
                    % checks to advance row counter and avoid overwriting
                    if all(~cellfun(@isempty,GABAmotion(GABAmotionRow,length(gaba_columns)+1:length(GABAmotion(1,:))))) % check to advance counter
                        GABAmotionRow = GABAmotionRow+1; % advance counter
                    elseif any(~cellfun(@isempty,GABAmotion(GABAmotionRow,length(gaba_columns)+1:length(GABAmotion(1,:))))) % check that data is not being overwritten
                        disp(['Avoided overwriting for asc file: ' data.experiment(expIdx).ascFilename{iRun} ' at row ' num2str(GABAmotionRow) ' of GABAmotion']) % message
                        GABAmotionRow = GABAmotionRow+1; % advance counter
                    end

                    % store metadata about eye tracking data
                    GABAmotion{GABAmotionRow,1} = gaba_subjects{iS}; % store subject name
                    GABAmotion{GABAmotionRow,2} = iRun; % store run number
                    GABAmotion{GABAmotionRow,3} = data.directory; % store directory
                    GABAmotion{GABAmotionRow,4} = data.experiment(expIdx).behavioralDataFilename{iRun}; % store behavioral data filename
                    GABAmotion{GABAmotionRow,5} = data.experiment(expIdx).ascFilename{iRun}; % store eyetracking data filename

                    % store stats by condition
                    for iCond = 1:length(motion_CondList)
                        if (size(data.experiment(expIdx).runStats{iRun},1) == 6 && (iCond == 5 || iCond == 6)) || ~istable(data.experiment(expIdx).runStats{iRun})
                            % conditions 5 and 6 (1 degree stimuli were
                            % not run for early GABA participants) 
                        else
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 6) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'noTrackTime')); % store no tracking time for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 7) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'fixationTime')); % store fixation time for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 8) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedFixationTime')); % store drift-corrected fixation time for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 9) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'saccadeTime')); % store saccade and blink time for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 10) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'distanceMean')); % store distance M for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 11) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'distanceSD')); % store distance SD for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 12) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'angleMean')); % store angle M for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 13) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'angleSD')); % store angle SD for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 14) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedDistanceMean')); % store drift-corrected distance M for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 15) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedDistanceSD')); % store drift-corrected distance SD for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 16) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedAngleMean')); % store drift-corrected angle M for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 17) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'driftCorrectedAngleSD')); % store drift-corrected angle SD for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 18) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeN')); % store number of saccades for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 19) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'bigSaccadeN')); % store number of big saccades for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 20) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'BlinkN')); % store number of saccades for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 21) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeAmplitudeMean')); % store saccade amplitude mean for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 22) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeAmplitudeSD')); % store saccade amplitude SD for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 23) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeVelocityMean')); % store saccade velocity mean for condition
                            GABAmotion(GABAmotionRow, length(statList) * (iCond - 1) + 24) = table2cell(data.experiment(expIdx).runStats{iRun}(motion_CondList{iCond},'SaccadeVelocitySD')); % store saccade velocity SD for condition
                        end
                    end
                end
            end
        end

        % get contrast data 
        expIdx = 2; % contrast
        if ~isempty(data.experiment(expIdx).ascFilename) && (iscell(data.experiment(expIdx).ascFilename) || ~isnan(data.experiment(expIdx).ascFilename)) % is there eye tracking data for this experiment?
            for iRun = 1:length(data.experiment(expIdx).runName)
                if onlyGood
                    if data.experiment(expIdx).quality{iRun} == 1
                        % checks to advance row counter and avoid overwriting
                        if all(~cellfun(@isempty,GABAcontrast(GABAcontrastRow,length(gaba_columns)+1:length(GABAcontrast(1,:))))) % check to advance counter
                            GABAcontrastRow = GABAcontrastRow+1; % advance counter
                        elseif any(~cellfun(@isempty,GABAcontrast(GABAcontrastRow,length(gaba_columns)+1:length(GABAcontrast(1,:))))) % check that data is not being overwritten
                            disp(['Avoided overwriting for asc file: ' data.experiment(expIdx).ascFilename{iRun} ' at row ' num2str(GABAcontrastRow) ' of GABAcontrast']) % message
                            GABAcontrastRow = GABAcontrastRow+1; % advance counter
                        end

                        % store metadata about eye tracking data
                        GABAcontrast{GABAcontrastRow,1} = gaba_subjects{iS}; % store subject name
                        GABAcontrast{GABAcontrastRow,2} = iRun; % store run number
                        GABAcontrast{GABAcontrastRow,3} = data.directory; % store directory
                        GABAcontrast{GABAcontrastRow,4} = data.experiment(expIdx).behavioralDataFilename{iRun}; % store behavioral data filename
                        GABAcontrast{GABAcontrastRow,5} = data.experiment(expIdx).ascFilename{iRun}; % store eyetracking data filename

                        % store stats by condition
                        for iCond = 1:length(contrast_gabaCondList)
                            if istable(data.experiment(expIdx).runStats{iRun})
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 6) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'noTrackTime')); % store no tracking time for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 7) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'fixationTime')); % store fixation time for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 8) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedFixationTime')); % store drift-corrected fixation time for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 9) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'saccadeTime')); % store saccade and blink time for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 10) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'distanceMean')); % store distance M for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 11) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'distanceSD')); % store distance SD for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 12) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'angleMean')); % store angle M for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 13) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'angleSD')); % store angle SD for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 14) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedDistanceMean')); % store drift-corrected distance M for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 15) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedDistanceSD')); % store drift-corrected distance SD for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 16) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedAngleMean')); % store drift-corrected angle M for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 17) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedAngleSD')); % store drift-corrected angle SD for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 18) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeN')); % store number of saccades for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 19) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'bigSaccadeN')); % store number of big saccades for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 20) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'BlinkN')); % store number of saccades for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 21) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeAmplitudeMean')); % store saccade amplitude mean for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 22) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeAmplitudeSD')); % store saccade amplitude SD for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 23) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeVelocityMean')); % store saccade velocity mean for condition
                                GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 24) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeVelocitySD')); % store saccade velocity SD for condition
                            end
                        end
                    end
                else
                    % checks to advance row counter and avoid overwriting
                    if all(~cellfun(@isempty,GABAcontrast(GABAcontrastRow,length(gaba_columns)+1:length(GABAcontrast(1,:))))) % check to advance counter
                        GABAcontrastRow = GABAcontrastRow+1; % advance counter
                    elseif any(~cellfun(@isempty,GABAcontrast(GABAcontrastRow,length(gaba_columns)+1:length(GABAcontrast(1,:))))) % check that data is not being overwritten
                        disp(['Avoided overwriting for asc file: ' data.experiment(expIdx).ascFilename{iRun} ' at row ' num2str(GABAcontrastRow) ' of GABAcontrast']) % message
                        GABAcontrastRow = GABAcontrastRow+1; % advance counter
                    end

                    % store metadata about eye tracking data
                    GABAcontrast{GABAcontrastRow,1} = gaba_subjects{iS}; % store subject name
                    GABAcontrast{GABAcontrastRow,2} = iRun; % store run number
                    GABAcontrast{GABAcontrastRow,3} = data.directory; % store directory
                    GABAcontrast{GABAcontrastRow,4} = data.experiment(expIdx).behavioralDataFilename{iRun}; % store behavioral data filename
                    GABAcontrast{GABAcontrastRow,5} = data.experiment(expIdx).ascFilename{iRun}; % store eyetracking data filename

                    % store stats by condition
                    for iCond = 1:length(contrast_gabaCondList)
                        if istable(data.experiment(expIdx).runStats{iRun})
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 6) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'noTrackTime')); % store no tracking time for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 7) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'fixationTime')); % store fixation time for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 8) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedFixationTime')); % store drift-corrected fixation time for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 9) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'saccadeTime')); % store saccade and blink time for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 10) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'distanceMean')); % store distance M for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 11) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'distanceSD')); % store distance SD for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 12) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'angleMean')); % store angle M for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 13) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'angleSD')); % store angle SD for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 14) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedDistanceMean')); % store drift-corrected distance M for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 15) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedDistanceSD')); % store drift-corrected distance SD for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 16) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedAngleMean')); % store drift-corrected angle M for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 17) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'driftCorrectedAngleSD')); % store drift-corrected angle SD for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 18) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeN')); % store number of saccades for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 19) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'bigSaccadeN')); % store number of big saccades for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 20) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'BlinkN')); % store number of saccades for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 21) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeAmplitudeMean')); % store saccade amplitude mean for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 22) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeAmplitudeSD')); % store saccade amplitude SD for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 23) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeVelocityMean')); % store saccade velocity mean for condition
                            GABAcontrast(GABAcontrastRow, length(statList) * (iCond - 1) + 24) = table2cell(data.experiment(expIdx).runStats{iRun}(contrast_gabaCondList{iCond},'SaccadeVelocitySD')); % store saccade velocity SD for condition
                        end
                    end
                end
            end
        end
    end
end

% build tables for Lorazepam data
for iS = 1:length(lz_subjects)
    % parse eye tracking data and save data structure for each subject individually
    clear eyetracking_data_file eyetracking_data_dir data
    eyetracking_data_file = [lz_subjects{iS} '_Psychophysics_EyeTracking_Data.mat'];
    eyetracking_data_dir = fullfile(top_dir,lz_dir,lz_subjects{iS});
    cd(eyetracking_data_dir);
    if ~exist(eyetracking_data_file,'file') % look for saved data
        data = AK_MotionContrastPsychophysics_EvaluateFixation(lz_subjects(iS));
        save(eyetracking_data_file,'data');
    elseif ~analyzeOnly
        disp(['loading ' eyetracking_data_file]) % message
        load(eyetracking_data_file,'data');
    end
    
    if ~analyzeOnly
        % get motion data 
        expIdx = 1; % motion
        if ~isempty(data.experiment(expIdx).ascFilename) && (iscell(data.experiment(expIdx).ascFilename) || ~isnan(data.experiment(expIdx).ascFilename)) % is there eye tracking data for this experiment?
            for iSess = 1:length(data.experiment(expIdx).runName(:,1));
                for iRun = 1:length(data.experiment(expIdx).runName(1,:))
                    if onlyGood
                        if data.experiment(expIdx).quality{iSess,iRun} == 1
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,LZmotion(LZmotionRow,length(gaba_columns)+1:length(LZmotion(1,:))))) % check to advance counter
                                LZmotionRow = LZmotionRow+1; % advance counter
                            elseif any(~cellfun(@isempty,LZmotion(LZmotionRow,length(gaba_columns)+1:length(LZmotion(1,:))))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.experiment(expIdx).ascFilename{iSess,iRun} ' at row ' num2str(LZmotionRow) ' of LZmotion']) % message
                                LZmotionRow = LZmotionRow+1; % advance counter
                            end

                            % store metadata about eye tracking data
                            LZmotion{LZmotionRow,1} = lz_subjects{iS}; % store subject name
                            LZmotion{LZmotionRow,2} = treatment{iS,iSess}; % store unblinded session id (drug vs placebo)
                            LZmotion{LZmotionRow,3} = iRun; % store run number
                            LZmotion{LZmotionRow,4} = data.directory{iSess}; % store directory
                            LZmotion{LZmotionRow,5} = data.experiment(expIdx).behavioralDataFilename{iSess,iRun}; % store behavioral data filename
                            LZmotion{LZmotionRow,6} = data.experiment(expIdx).ascFilename{iSess,iRun}; % store eyetracking data filename

                            % store stats by condition
                            for iCond = 1:length(motion_CondList)
                                if istable(data.experiment(expIdx).runStats{iSess,iRun})
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 7) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'noTrackTime')); % store no tracking time for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 8) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'fixationTime')); % store fixation time for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 9) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedFixationTime')); % store drift-corrected fixation time for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 10) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'saccadeTime')); % store saccade and blink time for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 11) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'distanceMean')); % store distance M for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 12) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'distanceSD')); % store distance SD for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 13) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'angleMean')); % store angle M for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 14) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'angleSD')); % store angle SD for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 15) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedDistanceMean')); % store drift-corrected distance M for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 16) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedDistanceSD')); % store drift-corrected distance SD for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 17) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedAngleMean')); % store drift-corrected angle M for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 18) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedAngleSD')); % store drift-corrected angle SD for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 19) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeN')); % store number of saccades for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 20) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'bigSaccadeN')); % store number of big saccades for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 21) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'BlinkN')); % store number of saccades for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 22) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeAmplitudeMean')); % store saccade amplitude mean for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 23) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeAmplitudeSD')); % store saccade amplitude SD for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 24) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeVelocityMean')); % store saccade velocity mean for condition
                                    LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 25) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeVelocitySD')); % store saccade velocity SD for condition
                                end
                            end
                        end
                    else
                        % checks to advance row counter and avoid overwriting
                        if all(~cellfun(@isempty,LZmotion(LZmotionRow,length(gaba_columns)+1:length(LZmotion(1,:))))) % check to advance counter
                            LZmotionRow = LZmotionRow+1; % advance counter
                        elseif any(~cellfun(@isempty,LZmotion(LZmotionRow,length(gaba_columns)+1:length(LZmotion(1,:))))) % check that data is not being overwritten
                            disp(['Avoided overwriting for asc file: ' data.experiment(expIdx).ascFilename{iSess,iRun} ' at row ' num2str(LZmotionRow) ' of LZmotion']) % message
                            LZmotionRow = LZmotionRow+1; % advance counter
                        end

                        % store metadata about eye tracking data
                        LZmotion{LZmotionRow,1} = lz_subjects{iS}; % store subject name
                        LZmotion{LZmotionRow,2} = treatment{iS,iSess}; % store unblinded session id (drug vs placebo)
                        LZmotion{LZmotionRow,3} = iRun; % store run number
                        LZmotion{LZmotionRow,4} = data.directory{iSess}; % store directory
                        LZmotion{LZmotionRow,5} = data.experiment(expIdx).behavioralDataFilename{iSess,iRun}; % store behavioral data filename
                        LZmotion{LZmotionRow,6} = data.experiment(expIdx).ascFilename{iSess,iRun}; % store eyetracking data filename

                        % store stats by condition
                        for iCond = 1:length(motion_CondList)
                            if istable(data.experiment(expIdx).runStats{iSess,iRun})
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 7) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'noTrackTime')); % store no tracking time for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 8) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'fixationTime')); % store fixation time for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 9) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedFixationTime')); % store drift-corrected fixation time for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 10) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'saccadeTime')); % store saccade and blink time for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 11) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'distanceMean')); % store distance M for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 12) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'distanceSD')); % store distance SD for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 13) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'angleMean')); % store angle M for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 14) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'angleSD')); % store angle SD for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 15) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedDistanceMean')); % store drift-corrected distance M for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 16) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedDistanceSD')); % store drift-corrected distance SD for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 17) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedAngleMean')); % store drift-corrected angle M for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 18) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'driftCorrectedAngleSD')); % store drift-corrected angle SD for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 19) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeN')); % store number of saccades for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 20) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'bigSaccadeN')); % store number of big saccades for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 21) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'BlinkN')); % store number of saccades for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 22) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeAmplitudeMean')); % store saccade amplitude mean for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 23) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeAmplitudeSD')); % store saccade amplitude SD for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 24) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeVelocityMean')); % store saccade velocity mean for condition
                                LZmotion(LZmotionRow, length(statList) * (iCond - 1) + 25) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(motion_CondList{iCond},'SaccadeVelocitySD')); % store saccade velocity SD for condition
                            end
                        end
                    end
                end
            end
        end

        % get contrast data 
        expIdx = 2; % contrast
        if ~isempty(data.experiment(expIdx).ascFilename) && (iscell(data.experiment(expIdx).ascFilename) || ~isnan(data.experiment(expIdx).ascFilename)) % is there eye tracking data for this experiment?
            for iSess = 1:length(data.experiment(expIdx).runName(:,1));
                for iRun = 1:length(data.experiment(expIdx).runName(1,:))
                    if onlyGood
                        if data.experiment(expIdx).quality{iSess,iRun} == 1
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,LZcontrast(LZcontrastRow,length(gaba_columns)+1:length(LZcontrast(1,:))))) % check to advance counter
                                LZcontrastRow = LZcontrastRow+1; % advance counter
                            elseif any(~cellfun(@isempty,LZcontrast(LZcontrastRow,length(gaba_columns)+1:length(LZcontrast(1,:))))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.experiment(expIdx).ascFilename{iSess,iRun} ' at row ' num2str(LZcontrastRow) ' of LZcontrast']) % message
                                LZcontrastRow = LZcontrastRow+1; % advance counter
                            end

                            % store metadata about eye tracking data
                            LZcontrast{LZcontrastRow,1} = lz_subjects{iS}; % store subject name
                            LZcontrast{LZcontrastRow,2} = treatment{iS,iSess}; % store unblinded session id (drug vs placebo)
                            LZcontrast{LZcontrastRow,3} = iRun; % store run number
                            LZcontrast{LZcontrastRow,4} = data.directory{iSess}; % store directory
                            LZcontrast{LZcontrastRow,5} = data.experiment(expIdx).behavioralDataFilename{iSess,iRun}; % store behavioral data filename
                            LZcontrast{LZcontrastRow,6} = data.experiment(expIdx).ascFilename{iSess,iRun}; % store eyetracking data filename

                            % store stats by condition
                            for iCond = 1:length(contrast_lzCondList)
                                if istable(data.experiment(expIdx).runStats{iSess,iRun})
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 7) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'noTrackTime')); % store no tracking time for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 8) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'fixationTime')); % store fixation time for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 9) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedFixationTime')); % store drift-corrected fixation time for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 10) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'saccadeTime')); % store saccade and blink time for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 11) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'distanceMean')); % store distance M for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 12) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'distanceSD')); % store distance SD for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 13) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'angleMean')); % store angle M for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 14) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'angleSD')); % store angle SD for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 15) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedDistanceMean')); % store drift-corrected distance M for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 16) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedDistanceSD')); % store drift-corrected distance SD for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 17) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedAngleMean')); % store drift-corrected angle M for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 18) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedAngleSD')); % store drift-corrected angle SD for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 19) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeN')); % store number of saccades for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 20) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'bigSaccadeN')); % store number of big saccades for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 21) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'BlinkN')); % store number of saccades for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 22) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeAmplitudeMean')); % store saccade amplitude mean for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 23) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeAmplitudeSD')); % store saccade amplitude SD for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 24) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeVelocityMean')); % store saccade velocity mean for condition
                                    LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 25) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeVelocitySD')); % store saccade velocity SD for condition
                                end
                            end
                        end
                    else
                        % checks to advance row counter and avoid overwriting
                        if all(~cellfun(@isempty,LZcontrast(LZcontrastRow,length(gaba_columns)+1:length(LZcontrast(1,:))))) % check to advance counter
                            LZcontrastRow = LZcontrastRow+1; % advance counter
                        elseif any(~cellfun(@isempty,LZcontrast(LZcontrastRow,length(gaba_columns)+1:length(LZcontrast(1,:))))) % check that data is not being overwritten
                            disp(['Avoided overwriting for asc file: ' data.experiment(expIdx).ascFilename{iSess,iRun} ' at row ' num2str(LZcontrastRow) ' of LZcontrast']) % message
                            LZcontrastRow = LZcontrastRow+1; % advance counter
                        end

                        % store metadata about eye tracking data
                        LZcontrast{LZcontrastRow,1} = lz_subjects{iS}; % store subject name
                        LZcontrast{LZcontrastRow,2} = treatment{iS,iSess}; % store unblinded session id (drug vs placebo)
                        LZcontrast{LZcontrastRow,3} = iRun; % store run number
                        LZcontrast{LZcontrastRow,4} = data.directory{iSess}; % store directory
                        LZcontrast{LZcontrastRow,5} = data.experiment(expIdx).behavioralDataFilename{iSess,iRun}; % store behavioral data filename
                        LZcontrast{LZcontrastRow,6} = data.experiment(expIdx).ascFilename{iSess,iRun}; % store eyetracking data filename

                        % store stats by condition
                        for iCond = 1:length(contrast_lzCondList)
                            if istable(data.experiment(expIdx).runStats{iSess,iRun})
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 7) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'noTrackTime')); % store no tracking time for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 8) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'fixationTime')); % store fixation time for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 9) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedFixationTime')); % store drift-corrected fixation time for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 10) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'saccadeTime')); % store saccade and blink time for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 11) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'distanceMean')); % store distance M for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 12) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'distanceSD')); % store distance SD for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 13) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'angleMean')); % store angle M for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 14) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'angleSD')); % store angle SD for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 15) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedDistanceMean')); % store drift-corrected distance M for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 16) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedDistanceSD')); % store drift-corrected distance SD for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 17) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedAngleMean')); % store drift-corrected angle M for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 18) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'driftCorrectedAngleSD')); % store drift-corrected angle SD for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 19) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeN')); % store number of saccades for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 20) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'bigSaccadeN')); % store number of big saccades for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 21) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'BlinkN')); % store number of saccades for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 22) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeAmplitudeMean')); % store saccade amplitude mean for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 23) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeAmplitudeSD')); % store saccade amplitude SD for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 24) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeVelocityMean')); % store saccade velocity mean for condition
                                LZcontrast(LZcontrastRow, length(statList) * (iCond - 1) + 25) = table2cell(data.experiment(expIdx).runStats{iSess,iRun}(contrast_lzCondList{iCond},'SaccadeVelocitySD')); % store saccade velocity SD for condition
                            end
                        end
                    end
                end
            end
        end
    end
end

if ~analyzeOnly
    % trim empty rows from end of each cell array
    GABAmotion(all(cellfun(@isempty,GABAmotion(:,1)),2),:) = [];
    GABAcontrast(all(cellfun(@isempty,GABAcontrast(:,1)),2),:) = [];
    LZmotion(all(cellfun(@isempty,LZmotion(:,1)),2),:) = [];
    LZcontrast(all(cellfun(@isempty,LZcontrast(:,1)),2),:) = [];

    %% write cell arrays to xls file spreadsheets

    success = zeros(4,1); % preallocate
    if onlyGood
        xlsFilename = fullfile(save_dir,'Psychophysics_EyeTracking_Data_Tables_OnlyGood.xlsx'); % name file for only subjects with good data
    else
        xlsFilename = fullfile(save_dir,'Psychophysics_EyeTracking_Data_Tables.xlsx'); % name file
    end

    success(1) = xlswrite(xlsFilename,GABAmotion,'GABA Motion');
    success(2) = xlswrite(xlsFilename,GABAcontrast,'GABA Contrast');
    success(3) = xlswrite(xlsFilename,LZmotion,'LZ Motion');
    success(4) = xlswrite(xlsFilename,LZcontrast,'LZ Contrast');
end
