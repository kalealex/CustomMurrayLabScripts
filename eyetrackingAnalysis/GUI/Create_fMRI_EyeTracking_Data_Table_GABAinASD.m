% Create_fMRI_EyeTracking_Data_Table_GABAinASD
% Murray_Lab_2016
% Created by AMK on 7/12/16, modefied to add block data 7/28/16, modefied
% to accomidate quality table 10/14/16, modefied to include ftap data
% 02/17/17

% switch: create table for good quality data only?
onlyGood = 1;

%% designate directory, subject, and logfile indentifying information

addpath(genpath('L:\MurrayLab\ASD\Data'));

top_dir = 'L:\MurrayLab\ASD\Data';
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab';
save_dir = 'L:\MurrayLab\DataTablesForGUI';
% strings for logfile/PRT file recognition and association
MTlocStr = {'MTlocalizer'};
V1locStr = {'V1localizer'};
V1_fixlocStr = {'V1localizer_fix'};
contrastStr = {'contrast'};
supStr = {'suppression'};
sumStr = {'summation'};
ftapStr = {'fingertapping'};

%% load eyetracking quality table to find subjects to include in analysis

load(fullfile(top_dir,'EyetrackingQuality.mat'));
subjects = qualityTable(2:end,1); % full list
subjects(cellfun(@isempty,subjects)) = []; % remove empty cells
subjects = unique(subjects); % remove repetitions from list

%% designate table structure and field names

% name columns of eye tracking data tables (condition name: statistic name)
summary_columns = {'subject','session','fMRI set#','logfile scan#'};
block_columns = [summary_columns(1:end) {'block#'}];
moving_columns = {'moving: no tracking time','moving: fix time','moving: drift corrected fix time','moving: saccade and blink time','moving: distance M','moving: distance SD','moving: angle M','moving: angle SD','moving: drift corrected distance M','moving: drift corrected distance SD','moving: drift corrected angle M','moving: drift corrected angle SD','moving: n saccades','moving: n big saccades','moving: n blinks','moving: saccade amplitude M','moving: saccade amplitude SD','moving: saccade velocity M','moving: saccade velocity SD'};
static_columns = {'static: no tracking time','static: fix time','static: drift corrected fix time','static: saccade and blink time','static: distance M','static: distance SD','static: angle M','static: angle SD','static: drift corrected distance M','static: drift corrected distance SD','static: drift corrected angle M','static: drift corrected angle SD','static: n saccades','static: n big saccades','static: n blinks','static: saccade amplitude M','static: saccade amplitude SD','static: saccade velocity M','static: saccade velocity SD'};
big_columns = {'big: no tracking time','big: fix time','big: drift corrected fix time','big: saccade and blink time','big: distance M','big: distance SD','big: angle M','big: angle SD','big: drift corrected distance M','big: drift corrected distance SD','big: drift corrected angle M','big: drift corrected angle SD','big: n saccades','big: n big saccades','big: n blinks','big: saccade amplitude M','big: saccade amplitude SD','big: saccade velocity M','big: saccade velocity SD'};
small_columns = {'small: no tracking time','small: fix time','small: drift corrected fix time','small: saccade and blink time','small: distance M','small: distance SD','small: angle M','small: angle SD','small: drift corrected distance M','small: drift corrected distance SD','small: drift corrected angle M','small: drift corrected angle SD','small: n saccades','small: n big saccades','small: n blinks','small: saccade amplitude M','small: saccade amplitude SD','small: saccade velocity M','small: saccade velocity SD'};
fix_columns = {'fix: no tracking time','fix: fix time','fix: drift corrected fix time','fix: saccade and blink time','fix: distance M','fix: distance SD','fix: angle M','fix: angle SD','fix: drift corrected distance M','fix: drift corrected distance SD','fix: drift corrected angle M','fix: drift corrected angle SD','fix: n saccades','fix: n big saccades','fix: n blinks','fix: saccade amplitude M','fix: saccade amplitude SD','fix: saccade velocity M','fix: saccade velocity SD'};
low_columns = {'low: no tracking time','low: fix time','low: drift corrected fix time','low: saccade and blink time','low: distance M','low: distance SD','low: angle M','low: angle SD','low: drift corrected distance M','low: drift corrected distance SD','low: drift corrected angle M','low: drift corrected angle SD','low: n saccades','low: n big saccades','low: n blinks','low: saccade amplitude M','low: saccade amplitude SD','low: saccade velocity M','low: saccade velocity SD'};
high_columns = {'high: no tracking time','high: fix time','high: drift corrected fix time','high: saccade and blink time','high: distance M','high: distance SD','high: angle M','high: angle SD','high: drift corrected distance M','high: drift corrected distance SD','high: drift corrected angle M','high: drift corrected angle SD','high: n saccades','high: n big saccades','high: n blinks','high: saccade amplitude M','high: saccade amplitude SD','high: saccade velocity M','high: saccade velocity SD'};
predictable_columns = {'predictable: no tracking time','predictable: fix time','predictable: drift corrected fix time','predictable: saccade and blink time','predictable: distance M','predictable: distance SD','predictable: angle M','predictable: angle SD','predictable: drift corrected distance M','predictable: drift corrected distance SD','predictable: drift corrected angle M','predictable: drift corrected angle SD','predictable: n saccades','predictable: n big saccades','predictable: n blinks','predictable: saccade amplitude M','predictable: saccade amplitude SD','predictable: saccade velocity M','predictable: saccade velocity SD'};
variable_columns = {'variable: no tracking time','variable: fix time','variable: drift corrected fix time','variable: saccade and blink time','variable: distance M','variable: distance SD','variable: angle M','variable: angle SD','variable: drift corrected distance M','variable: drift corrected distance SD','variable: drift corrected angle M','variable: drift corrected angle SD','variable: n saccades','variable: n big saccades','variable: n blinks','variable: saccade amplitude M','variable: saccade amplitude SD','variable: saccade velocity M','variable: saccade velocity SD'};
all_columns = {'all: no tracking time','all: fix time','all: drift corrected fix time','all: saccade and blink time','all: distance M','all: distance SD','all: angle M','all: angle SD','all: drift corrected distance M','all: drift corrected distance SD','all: drift corrected angle M','all: drift corrected angle SD','all: n saccades','all: n big saccades','all: n blinks','all: saccade amplitude M','all: saccade amplitude SD','all: saccade velocity M','all: saccade velocity SD'};

% create cell arrays to store data and save into spreadsheets
MTloc = cell(3*length(subjects)+50,length(summary_columns)+length(moving_columns)+length(static_columns)+length(all_columns)); % 50 rows for buffer
V1loc = cell(3*length(subjects)+50,length(summary_columns)+length(big_columns)+length(small_columns)+length(all_columns)); % 50 rows for buffer
V1_fixloc = cell(3*length(subjects)+50,length(summary_columns)+length(fix_columns)+length(small_columns)+length(all_columns)); % 50 rows for buffer
contrast = cell(6*length(subjects)+100,length(summary_columns)+length(fix_columns)+length(low_columns)+length(high_columns)+length(all_columns)); % 100 rows for buffer
sup = cell(6*length(subjects)+100,length(summary_columns)+length(small_columns)+length(big_columns)+length(all_columns)); % 100 rows for buffer; sup and sum have the same condition names but separate cell arrays
sum = cell(6*length(subjects)+100,length(summary_columns)+length(small_columns)+length(big_columns)+length(all_columns)); % 100 rows for buffer
ftap = cell(6*length(subjects)+150,length(summary_columns)+length(fix_columns)+length(predictable_columns)+length(variable_columns)+length(all_columns)); % 150 rows for buffer

bl_MTloc = cell(3*length(subjects)+50,length(block_columns)+length(moving_columns)+length(static_columns)+length(all_columns)); % 50 rows for buffer
bl_V1loc = cell(3*length(subjects)+50,length(block_columns)+length(big_columns)+length(small_columns)+length(all_columns)); % 50 rows for buffer
bl_V1_fixloc = cell(3*length(subjects)+50,length(block_columns)+length(fix_columns)+length(small_columns)+length(all_columns)); % 50 rows for buffer
bl_contrast = cell(6*length(subjects)+100,length(block_columns)+length(fix_columns)+length(low_columns)+length(high_columns)+length(all_columns)); % 100 rows for buffer
bl_sup = cell(6*length(subjects)+100,length(block_columns)+length(small_columns)+length(big_columns)+length(all_columns)); % 100 rows for buffer; sup and sum have the same condition names but separate cell arrays
bl_sum = cell(6*length(subjects)+100,length(block_columns)+length(small_columns)+length(big_columns)+length(all_columns)); % 100 rows for buffer
bl_ftap = cell(6*length(subjects)+150,length(block_columns)+length(fix_columns)+length(predictable_columns)+length(variable_columns)+length(all_columns)); % 150 rows for buffer

% add column names to cell arrays
MTloc(1,:) = [summary_columns,moving_columns,static_columns,all_columns];
V1loc(1,:) = [summary_columns,big_columns,small_columns,all_columns];
V1_fixloc(1,:) = [summary_columns,fix_columns,small_columns,all_columns];
contrast(1,:) = [summary_columns,fix_columns,low_columns,high_columns,all_columns];
sup(1,:) = [summary_columns,small_columns,big_columns,all_columns];
sum(1,:) = [summary_columns,small_columns,big_columns,all_columns];
ftap(1,:) = [summary_columns,fix_columns,predictable_columns,variable_columns,all_columns];

bl_MTloc(1,:) = [block_columns,moving_columns,static_columns,all_columns];
bl_V1loc(1,:) = [block_columns,big_columns,small_columns,all_columns];
bl_V1_fixloc(1,:) = [block_columns,fix_columns,small_columns,all_columns];
bl_contrast(1,:) = [block_columns,fix_columns,low_columns,high_columns,all_columns];
bl_sup(1,:) = [block_columns,small_columns,big_columns,all_columns];
bl_sum(1,:) = [block_columns,small_columns,big_columns,all_columns];
bl_ftap(1,:) = [block_columns,fix_columns,predictable_columns,variable_columns,all_columns];

% create counters for rows of each cell array
MTlocRow = 2;
V1locRow = 2;
V1_fixlocRow = 2;
contrastRow = 2;
supRow = 2;
sumRow = 2;
ftapRow = 2;

bl_MTlocRow = 2;
bl_V1locRow = 2;
bl_V1_fixlocRow = 2;
bl_contrastRow = 2;
bl_supRow = 2;
bl_sumRow = 2;
bl_ftapRow = 2;

%% pull data from data structure into appropriate spreadsheets

for iS = 1:length(subjects) % cycle through subjects 
    % parse eye tracking data and save data structure for each subject individually
    clear eyetracking_data_filename eyetracking_data_dir data
    eyetracking_data_filename = [subjects{iS} '_fMRI_EyeTracking_Data.mat'];
    eyetracking_data_dir = fullfile(top_dir,subjects{iS},eyetracking_data_filename);
    if ~exist(eyetracking_data_dir,'file') % look for saved data
        data = AK_GABAinASD_fMRI_EvaluateFixation(subjects(iS));
        save(eyetracking_data_dir,'data');
    else
        disp(['loading ' eyetracking_data_dir]) % message
        load(eyetracking_data_dir,'data');
    end
    
    if isfield(data,'ascFilename') % will exclude disqualified subjects from table
        iF = 1; % reset fMRIsession counter
        while iF <= length(data.ascFilename(:,1)) % cycle through fMRI sessions
            for iA = 1:length(data.ascFilename(iF,:)) % cycle through asc files
                clear scanNfindStr scanNindex scanN
                scanNfindStr = [subjects{iS} '_']; % string to be used to find scan# and logfile name
                scanNindex = strfind(data.ascFilename{iF,iA},scanNfindStr);               
                % find and store scan# as variable scanN
                if ~isnan(str2double(data.ascFilename{iF,iA}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1))) % check for two digit scan#
                    scanN = str2double(data.ascFilename{iF,iA}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1));
                elseif ~isnan(str2double(data.ascFilename{iF,iA}(scanNindex+length(scanNfindStr)))) % check for one digit scan#
                    scanN = str2double(data.ascFilename{iF,iA}(scanNindex+length(scanNfindStr)));
                end
                % either look at only the good quality eye tracking data or
                % all eye tracking data that is actually data (not just
                % noise) based on qualityTable
                switch onlyGood
                    case 0
                        % determine which type of scan iA corresponds to
                        if AK_whichPattern(data.logFilename{iF,iA},MTlocStr,true)==1 && istable(data.summaryStats{iF,iA}) % is this an MT localizer?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,MTloc(MTlocRow,5:42))) % check to advance counter
                                MTlocRow = MTlocRow+1; % advance counter
                            elseif any(~cellfun(@isempty,MTloc(MTlocRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(MTlocRow) ' of MTloc']) % message
                                MTlocRow = MTlocRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            MTloc{MTlocRow,1} = subjects{iS}; % store subject name
                            MTloc{MTlocRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            MTloc{MTlocRow,3} = data.setN{iF,iA}; % store set number
                            MTloc{MTlocRow,4} = num2str(scanN); % store scan#

                            MTloc(MTlocRow,5) = table2cell(data.summaryStats{iF,iA}('moving','noTrackTime')); % store no tracking time for moving condition
                            MTloc(MTlocRow,6) = table2cell(data.summaryStats{iF,iA}('moving','fixationTime')); % store fixation time for moving condition
                            MTloc(MTlocRow,7) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedFixationTime')); % store drift-corrected fixation time for moving condition
                            MTloc(MTlocRow,8) = table2cell(data.summaryStats{iF,iA}('moving','saccadeTime')); % store saccade and blink time for moving condition
                            MTloc(MTlocRow,9) = table2cell(data.summaryStats{iF,iA}('moving','distanceMean')); % store distance M for moving condition
                            MTloc(MTlocRow,10) = table2cell(data.summaryStats{iF,iA}('moving','distanceSD')); % store distance SD for moving condition
                            MTloc(MTlocRow,11) = table2cell(data.summaryStats{iF,iA}('moving','angleMean')); % store angle M for moving condition
                            MTloc(MTlocRow,12) = table2cell(data.summaryStats{iF,iA}('moving','angleSD')); % store angle SD for moving condition
                            MTloc(MTlocRow,13) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedDistanceMean')); % store drift-corrected distance M for moving condition
                            MTloc(MTlocRow,14) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedDistanceSD')); % store drift-corrected distance SD for moving condition
                            MTloc(MTlocRow,15) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedAngleMean')); % store drift-corrected angle M for moving condition
                            MTloc(MTlocRow,16) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedAngleSD')); % store drift-corrected angle SD for moving condition
                            MTloc(MTlocRow,17) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeN')); % store number of saccades for moving condition
                            MTloc(MTlocRow,18) = table2cell(data.summaryStats{iF,iA}('moving','bigSaccadeN')); % store number of big saccades for moving condition
                            MTloc(MTlocRow,19) = table2cell(data.summaryStats{iF,iA}('moving','BlinkN')); % store number of saccades for moving condition
                            MTloc(MTlocRow,20) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeAmplitudeMean')); % store saccade amplitude mean for moving condition
                            MTloc(MTlocRow,21) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeAmplitudeSD')); % store saccade amplitude SD for moving condition
                            MTloc(MTlocRow,22) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeVelocityMean')); % store saccade velocity mean for moving condition
                            MTloc(MTlocRow,23) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeVelocitySD')); % store saccade velocity SD for moving condition

                            MTloc(MTlocRow,24) = table2cell(data.summaryStats{iF,iA}('static','noTrackTime')); % store no tracking time for static condition
                            MTloc(MTlocRow,25) = table2cell(data.summaryStats{iF,iA}('static','fixationTime')); % store fixation time for static condition
                            MTloc(MTlocRow,26) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedFixationTime')); % store drift-corrected fixation time for static condition
                            MTloc(MTlocRow,27) = table2cell(data.summaryStats{iF,iA}('static','saccadeTime')); % store saccade and blink time for static condition
                            MTloc(MTlocRow,28) = table2cell(data.summaryStats{iF,iA}('static','distanceMean')); % store distance M for static condition
                            MTloc(MTlocRow,29) = table2cell(data.summaryStats{iF,iA}('static','distanceSD')); % store distance SD for static condition
                            MTloc(MTlocRow,30) = table2cell(data.summaryStats{iF,iA}('static','angleMean')); % store angle M for static condition
                            MTloc(MTlocRow,31) = table2cell(data.summaryStats{iF,iA}('static','angleSD')); % store angle SD for static condition
                            MTloc(MTlocRow,32) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedDistanceMean')); % store drift-corrected distance M for static condition
                            MTloc(MTlocRow,33) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedDistanceSD')); % store drift-corrected distance SD for static condition
                            MTloc(MTlocRow,34) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedAngleMean')); % store drift-corrected angle M for static condition
                            MTloc(MTlocRow,35) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedAngleSD')); % store drift-corrected angle SD for static condition
                            MTloc(MTlocRow,36) = table2cell(data.summaryStats{iF,iA}('static','SaccadeN')); % store number of saccades for static condition
                            MTloc(MTlocRow,37) = table2cell(data.summaryStats{iF,iA}('static','bigSaccadeN')); % store number of big saccades for static condition
                            MTloc(MTlocRow,38) = table2cell(data.summaryStats{iF,iA}('static','BlinkN')); % store number of saccades for static condition
                            MTloc(MTlocRow,39) = table2cell(data.summaryStats{iF,iA}('static','SaccadeAmplitudeMean')); % store saccade amplitude mean for static condition
                            MTloc(MTlocRow,40) = table2cell(data.summaryStats{iF,iA}('static','SaccadeAmplitudeSD')); % store saccade amplitude SD for static condition
                            MTloc(MTlocRow,41) = table2cell(data.summaryStats{iF,iA}('static','SaccadeVelocityMean')); % store saccade velocity mean for static condition
                            MTloc(MTlocRow,42) = table2cell(data.summaryStats{iF,iA}('static','SaccadeVelocitySD')); % store saccade velocity SD for static condition
                            
                            MTloc(MTlocRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            MTloc(MTlocRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            MTloc(MTlocRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            MTloc(MTlocRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            MTloc(MTlocRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            MTloc(MTlocRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            MTloc(MTlocRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            MTloc(MTlocRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            MTloc(MTlocRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            MTloc(MTlocRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            MTloc(MTlocRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            MTloc(MTlocRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            MTloc(MTlocRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            MTloc(MTlocRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            MTloc(MTlocRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            MTloc(MTlocRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            MTloc(MTlocRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            MTloc(MTlocRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            MTloc(MTlocRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_MTloc{bl_MTlocRow,1} = subjects{iS}; % store subject name
                                bl_MTloc{bl_MTlocRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_MTloc{bl_MTlocRow,3} = data.setN{iF,iA}; % store set number
                                bl_MTloc{bl_MTlocRow,4} = num2str(scanN); % store scan#
                                bl_MTloc(bl_MTlocRow,5) = {iB}; % store block#

                                bl_MTloc(bl_MTlocRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for moving condition
                                bl_MTloc(bl_MTlocRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for moving condition
                                bl_MTloc(bl_MTlocRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for moving condition
                                bl_MTloc(bl_MTlocRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for moving condition
                                bl_MTloc(bl_MTlocRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''distanceMean''}{1}(iB)}','cell'); % store distance M for moving condition
                                bl_MTloc(bl_MTlocRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for moving condition
                                bl_MTloc(bl_MTlocRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''angleMean''}{1}(iB)}','cell'); % store angle M for moving condition
                                bl_MTloc(bl_MTlocRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''angleSD''}{1}(iB)}','cell'); % store angle SD for moving condition
                                bl_MTloc(bl_MTlocRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for moving condition
                                bl_MTloc(bl_MTlocRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for moving condition
                                bl_MTloc(bl_MTlocRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for moving condition
                                bl_MTloc(bl_MTlocRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for moving condition
                                bl_MTloc(bl_MTlocRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for moving condition
                                bl_MTloc(bl_MTlocRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for moving condition
                                bl_MTloc(bl_MTlocRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for moving condition
                                bl_MTloc(bl_MTlocRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for moving condition
                                bl_MTloc(bl_MTlocRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for moving condition
                                bl_MTloc(bl_MTlocRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for moving condition
                                bl_MTloc(bl_MTlocRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for moving condition

                                bl_MTloc(bl_MTlocRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for static condition
                                bl_MTloc(bl_MTlocRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for static condition
                                bl_MTloc(bl_MTlocRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for static condition
                                bl_MTloc(bl_MTlocRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for static condition
                                bl_MTloc(bl_MTlocRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''distanceMean''}{1}(iB)}','cell'); % store distance M for static condition
                                bl_MTloc(bl_MTlocRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for static condition
                                bl_MTloc(bl_MTlocRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''angleMean''}{1}(iB)}','cell'); % store angle M for static condition
                                bl_MTloc(bl_MTlocRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''angleSD''}{1}(iB)}','cell'); % store angle SD for static condition
                                bl_MTloc(bl_MTlocRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for static condition
                                bl_MTloc(bl_MTlocRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for static condition
                                bl_MTloc(bl_MTlocRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for static condition
                                bl_MTloc(bl_MTlocRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for static condition
                                bl_MTloc(bl_MTlocRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for static condition
                                bl_MTloc(bl_MTlocRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for static condition
                                bl_MTloc(bl_MTlocRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for static condition
                                bl_MTloc(bl_MTlocRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for static condition
                                bl_MTloc(bl_MTlocRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for static condition
                                bl_MTloc(bl_MTlocRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for static condition
                                bl_MTloc(bl_MTlocRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for static condition

                                bl_MTlocRow = bl_MTlocRow+1; % advance counter
                            end

                        elseif AK_whichPattern(data.logFilename{iF,iA},V1_fixlocStr,true)==1 && istable(data.summaryStats{iF,iA}) % is this an V1_fix localizer?; must check before V1 localizer because of overlap in strings
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,V1_fixloc(V1_fixlocRow,5:42))) % check to advance counter
                                V1_fixlocRow = V1_fixlocRow+1; % advance counter
                            elseif any(~cellfun(@isempty,V1_fixloc(V1_fixlocRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(V1_fixlocRow) ' of V1_fixloc']) % message
                                V1_fixlocRow = V1_fixlocRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            V1_fixloc{V1_fixlocRow,1} = subjects{iS}; % store subject name
                            V1_fixloc{V1_fixlocRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            V1_fixloc{V1_fixlocRow,3} = data.setN{iF,iA}; % store set number
                            V1_fixloc{V1_fixlocRow,4} = num2str(scanN); % store scan#

                            V1_fixloc(V1_fixlocRow,5) = table2cell(data.summaryStats{iF,iA}('surr','noTrackTime')); % store no tracking time for fix condition
                            V1_fixloc(V1_fixlocRow,6) = table2cell(data.summaryStats{iF,iA}('surr','fixationTime')); % store fixation time for fix condition
                            V1_fixloc(V1_fixlocRow,7) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedFixationTime')); % store drift-corrected fixation time for fix condition
                            V1_fixloc(V1_fixlocRow,8) = table2cell(data.summaryStats{iF,iA}('surr','saccadeTime')); % store saccade and blink time for fix condition
                            V1_fixloc(V1_fixlocRow,9) = table2cell(data.summaryStats{iF,iA}('surr','distanceMean')); % store distance M for fix condition
                            V1_fixloc(V1_fixlocRow,10) = table2cell(data.summaryStats{iF,iA}('surr','distanceSD')); % store distance SD for fix condition
                            V1_fixloc(V1_fixlocRow,11) = table2cell(data.summaryStats{iF,iA}('surr','angleMean')); % store angle M for fix condition
                            V1_fixloc(V1_fixlocRow,12) = table2cell(data.summaryStats{iF,iA}('surr','angleSD')); % store angle SD for fix condition
                            V1_fixloc(V1_fixlocRow,13) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedDistanceMean')); % store drift-corrected distance M for fix condition
                            V1_fixloc(V1_fixlocRow,14) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedDistanceSD')); % store drift-corrected distance SD for fix condition
                            V1_fixloc(V1_fixlocRow,15) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedAngleMean')); % store drift-corrected angle M for fix condition
                            V1_fixloc(V1_fixlocRow,16) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedAngleSD')); % store drift-corrected angle SD for fix condition
                            V1_fixloc(V1_fixlocRow,17) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeN')); % store number of saccades for fix condition
                            V1_fixloc(V1_fixlocRow,18) = table2cell(data.summaryStats{iF,iA}('surr','bigSaccadeN')); % store number of big saccades for fix condition
                            V1_fixloc(V1_fixlocRow,19) = table2cell(data.summaryStats{iF,iA}('surr','BlinkN')); % store number of saccades for fix condition
                            V1_fixloc(V1_fixlocRow,20) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeAmplitudeMean')); % store saccade amplitude mean for fix condition
                            V1_fixloc(V1_fixlocRow,21) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeAmplitudeSD')); % store saccade amplitude SD for fix condition
                            V1_fixloc(V1_fixlocRow,22) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeVelocityMean')); % store saccade velocity mean for fix condition
                            V1_fixloc(V1_fixlocRow,23) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeVelocitySD')); % store saccade velocity SD for fix condition

                            V1_fixloc(V1_fixlocRow,24) = table2cell(data.summaryStats{iF,iA}('tar','noTrackTime')); % store no tracking time for small condition
                            V1_fixloc(V1_fixlocRow,25) = table2cell(data.summaryStats{iF,iA}('tar','fixationTime')); % store fixation time for small condition
                            V1_fixloc(V1_fixlocRow,26) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedFixationTime')); % store drift-corrected fixation time for small condition
                            V1_fixloc(V1_fixlocRow,27) = table2cell(data.summaryStats{iF,iA}('tar','saccadeTime')); % store saccade and blink time for small condition
                            V1_fixloc(V1_fixlocRow,28) = table2cell(data.summaryStats{iF,iA}('tar','distanceMean')); % store distance M for small condition
                            V1_fixloc(V1_fixlocRow,29) = table2cell(data.summaryStats{iF,iA}('tar','distanceSD')); % store distance SD for small condition
                            V1_fixloc(V1_fixlocRow,30) = table2cell(data.summaryStats{iF,iA}('tar','angleMean')); % store angle M for small condition
                            V1_fixloc(V1_fixlocRow,31) = table2cell(data.summaryStats{iF,iA}('tar','angleSD')); % store angle SD for small condition
                            V1_fixloc(V1_fixlocRow,32) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedDistanceMean')); % store drift-corrected distance M for small condition
                            V1_fixloc(V1_fixlocRow,33) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedDistanceSD')); % store drift-corrected distance SD for small condition
                            V1_fixloc(V1_fixlocRow,34) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedAngleMean')); % store drift-corrected angle M for small condition
                            V1_fixloc(V1_fixlocRow,35) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedAngleSD')); % store drift-corrected angle SD for small condition
                            V1_fixloc(V1_fixlocRow,36) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeN')); % store number of saccades for small condition
                            V1_fixloc(V1_fixlocRow,37) = table2cell(data.summaryStats{iF,iA}('tar','bigSaccadeN')); % store number of big saccades for small condition
                            V1_fixloc(V1_fixlocRow,38) = table2cell(data.summaryStats{iF,iA}('tar','BlinkN')); % store number of saccades for small condition
                            V1_fixloc(V1_fixlocRow,39) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeAmplitudeMean')); % store saccade amplitude mean for small condition
                            V1_fixloc(V1_fixlocRow,40) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeAmplitudeSD')); % store saccade amplitude SD for small condition
                            V1_fixloc(V1_fixlocRow,41) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeVelocityMean')); % store saccade velocity mean for small condition
                            V1_fixloc(V1_fixlocRow,42) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeVelocitySD')); % store saccade velocity SD for small condition

                            V1_fixloc(V1_fixlocRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            V1_fixloc(V1_fixlocRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            V1_fixloc(V1_fixlocRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            V1_fixloc(V1_fixlocRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            V1_fixloc(V1_fixlocRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            V1_fixloc(V1_fixlocRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            V1_fixloc(V1_fixlocRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            V1_fixloc(V1_fixlocRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            V1_fixloc(V1_fixlocRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            V1_fixloc(V1_fixlocRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            V1_fixloc(V1_fixlocRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            V1_fixloc(V1_fixlocRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            V1_fixloc(V1_fixlocRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            V1_fixloc(V1_fixlocRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            V1_fixloc(V1_fixlocRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            V1_fixloc(V1_fixlocRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            V1_fixloc(V1_fixlocRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            V1_fixloc(V1_fixlocRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            V1_fixloc(V1_fixlocRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_V1_fixloc{bl_V1_fixlocRow,1} = subjects{iS}; % store subject name
                                bl_V1_fixloc{bl_V1_fixlocRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_V1_fixloc{bl_V1_fixlocRow,3} = data.setN{iF,iA}; % store set number
                                bl_V1_fixloc{bl_V1_fixlocRow,4} = num2str(scanN); % store scan#
                                bl_V1_fixloc(bl_V1_fixlocRow,5) = {iB}; % store block#

                                bl_V1_fixloc(bl_V1_fixlocRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''distanceMean''}{1}(iB)}','cell'); % store distance M for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''angleMean''}{1}(iB)}','cell'); % store angle M for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''angleSD''}{1}(iB)}','cell'); % store angle SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for fix condition

                                bl_V1_fixloc(bl_V1_fixlocRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''distanceMean''}{1}(iB)}','cell'); % store distance M for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''angleMean''}{1}(iB)}','cell'); % store angle M for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''angleSD''}{1}(iB)}','cell'); % store angle SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for small condition

                                bl_V1_fixlocRow = bl_V1_fixlocRow+1; % advance counter
                            end

                        elseif AK_whichPattern(data.logFilename{iF,iA},V1locStr,true)==1 && istable(data.summaryStats{iF,iA}) % is this an V1 localizer?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,V1loc(V1locRow,5:42))) % check to advance counter
                                V1locRow = V1locRow+1; % advance counter
                            elseif any(~cellfun(@isempty,V1loc(V1locRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(V1locRow) ' of V1loc']) % message
                                V1locRow = V1locRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            V1loc{V1locRow,1} = subjects{iS}; % store subject name
                            V1loc{V1locRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            V1loc{V1locRow,3} = data.setN{iF,iA}; % store set number
                            V1loc{V1locRow,4} = num2str(scanN); % store scan#

                            V1loc(V1locRow,5) = table2cell(data.summaryStats{iF,iA}('surr','noTrackTime')); % store no tracking time for big condition
                            V1loc(V1locRow,6) = table2cell(data.summaryStats{iF,iA}('surr','fixationTime')); % store fixation time for big condition
                            V1loc(V1locRow,7) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedFixationTime')); % store drift-corrected fixation time for big condition
                            V1loc(V1locRow,8) = table2cell(data.summaryStats{iF,iA}('surr','saccadeTime')); % store saccade and blink time for big condition
                            V1loc(V1locRow,9) = table2cell(data.summaryStats{iF,iA}('surr','distanceMean')); % store distance M for big condition
                            V1loc(V1locRow,10) = table2cell(data.summaryStats{iF,iA}('surr','distanceSD')); % store distance SD for big condition
                            V1loc(V1locRow,11) = table2cell(data.summaryStats{iF,iA}('surr','angleMean')); % store angle M for big condition
                            V1loc(V1locRow,12) = table2cell(data.summaryStats{iF,iA}('surr','angleSD')); % store angle SD for big condition
                            V1loc(V1locRow,13) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedDistanceMean')); % store drift-corrected distance M for big condition
                            V1loc(V1locRow,14) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedDistanceSD')); % store drift-corrected distance SD for big condition
                            V1loc(V1locRow,15) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedAngleMean')); % store drift-corrected angle M for big condition
                            V1loc(V1locRow,16) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedAngleSD')); % store drift-corrected angle SD for big condition
                            V1loc(V1locRow,17) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeN')); % store number of saccades for big condition
                            V1loc(V1locRow,18) = table2cell(data.summaryStats{iF,iA}('surr','bigSaccadeN')); % store number of big saccades for big condition
                            V1loc(V1locRow,19) = table2cell(data.summaryStats{iF,iA}('surr','BlinkN')); % store number of saccades for big condition
                            V1loc(V1locRow,20) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeAmplitudeMean')); % store saccade amplitude mean for big condition
                            V1loc(V1locRow,21) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeAmplitudeSD')); % store saccade amplitude SD for big condition
                            V1loc(V1locRow,22) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeVelocityMean')); % store saccade velocity mean for big condition
                            V1loc(V1locRow,23) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeVelocitySD')); % store saccade velocity SD for big condition

                            V1loc(V1locRow,24) = table2cell(data.summaryStats{iF,iA}('tar','noTrackTime')); % store no tracking time for small condition
                            V1loc(V1locRow,25) = table2cell(data.summaryStats{iF,iA}('tar','fixationTime')); % store fixation time for small condition
                            V1loc(V1locRow,26) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedFixationTime')); % store drift-corrected fixation time for small condition
                            V1loc(V1locRow,27) = table2cell(data.summaryStats{iF,iA}('tar','saccadeTime')); % store saccade and blink time for small condition
                            V1loc(V1locRow,28) = table2cell(data.summaryStats{iF,iA}('tar','distanceMean')); % store distance M for small condition
                            V1loc(V1locRow,29) = table2cell(data.summaryStats{iF,iA}('tar','distanceSD')); % store distance SD for small condition
                            V1loc(V1locRow,30) = table2cell(data.summaryStats{iF,iA}('tar','angleMean')); % store angle M for small condition
                            V1loc(V1locRow,31) = table2cell(data.summaryStats{iF,iA}('tar','angleSD')); % store angle SD for small condition
                            V1loc(V1locRow,32) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedDistanceMean')); % store drift-corrected distance M for small condition
                            V1loc(V1locRow,33) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedDistanceSD')); % store drift-corrected distance SD for small condition
                            V1loc(V1locRow,34) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedAngleMean')); % store drift-corrected angle M for small condition
                            V1loc(V1locRow,35) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedAngleSD')); % store drift-corrected angle SD for small condition
                            V1loc(V1locRow,36) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeN')); % store number of saccades for small condition
                            V1loc(V1locRow,37) = table2cell(data.summaryStats{iF,iA}('tar','bigSaccadeN')); % store number of big saccades for small condition
                            V1loc(V1locRow,38) = table2cell(data.summaryStats{iF,iA}('tar','BlinkN')); % store number of saccades for small condition
                            V1loc(V1locRow,39) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeAmplitudeMean')); % store saccade amplitude mean for small condition
                            V1loc(V1locRow,40) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeAmplitudeSD')); % store saccade amplitude SD for small condition
                            V1loc(V1locRow,41) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeVelocityMean')); % store saccade velocity mean for small condition
                            V1loc(V1locRow,42) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeVelocitySD')); % store saccade velocity SD for small condition

                            V1loc(V1locRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            V1loc(V1locRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            V1loc(V1locRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            V1loc(V1locRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            V1loc(V1locRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            V1loc(V1locRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            V1loc(V1locRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            V1loc(V1locRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            V1loc(V1locRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            V1loc(V1locRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            V1loc(V1locRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            V1loc(V1locRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            V1loc(V1locRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            V1loc(V1locRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            V1loc(V1locRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            V1loc(V1locRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            V1loc(V1locRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            V1loc(V1locRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            V1loc(V1locRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_V1loc{bl_V1locRow,1} = subjects{iS}; % store subject name
                                bl_V1loc{bl_V1locRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_V1loc{bl_V1locRow,3} = data.setN{iF,iA}; % store set number
                                bl_V1loc{bl_V1locRow,4} = num2str(scanN); % store scan#
                                bl_V1loc(bl_V1locRow,5) = {iB}; % store block#

                                bl_V1loc(bl_V1locRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for big condition
                                bl_V1loc(bl_V1locRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for big condition
                                bl_V1loc(bl_V1locRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for big condition
                                bl_V1loc(bl_V1locRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for big condition
                                bl_V1loc(bl_V1locRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''distanceMean''}{1}(iB)}','cell'); % store distance M for big condition
                                bl_V1loc(bl_V1locRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for big condition
                                bl_V1loc(bl_V1locRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''angleMean''}{1}(iB)}','cell'); % store angle M for big condition
                                bl_V1loc(bl_V1locRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''angleSD''}{1}(iB)}','cell'); % store angle SD for big condition
                                bl_V1loc(bl_V1locRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for big condition
                                bl_V1loc(bl_V1locRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for big condition
                                bl_V1loc(bl_V1locRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for big condition
                                bl_V1loc(bl_V1locRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for big condition
                                bl_V1loc(bl_V1locRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for big condition
                                bl_V1loc(bl_V1locRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for big condition
                                bl_V1loc(bl_V1locRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for big condition
                                bl_V1loc(bl_V1locRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for big condition
                                bl_V1loc(bl_V1locRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for big condition
                                bl_V1loc(bl_V1locRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for big condition
                                bl_V1loc(bl_V1locRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for big condition

                                bl_V1loc(bl_V1locRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for small condition
                                bl_V1loc(bl_V1locRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for small condition
                                bl_V1loc(bl_V1locRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for small condition
                                bl_V1loc(bl_V1locRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for small condition
                                bl_V1loc(bl_V1locRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''distanceMean''}{1}(iB)}','cell'); % store distance M for small condition
                                bl_V1loc(bl_V1locRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for small condition
                                bl_V1loc(bl_V1locRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''angleMean''}{1}(iB)}','cell'); % store angle M for small condition
                                bl_V1loc(bl_V1locRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''angleSD''}{1}(iB)}','cell'); % store angle SD for small condition
                                bl_V1loc(bl_V1locRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for small condition
                                bl_V1loc(bl_V1locRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for small condition
                                bl_V1loc(bl_V1locRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for small condition
                                bl_V1loc(bl_V1locRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for small condition
                                bl_V1loc(bl_V1locRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for small condition
                                bl_V1loc(bl_V1locRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for small condition
                                bl_V1loc(bl_V1locRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for small condition
                                bl_V1loc(bl_V1locRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for small condition
                                bl_V1loc(bl_V1locRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for small condition
                                bl_V1loc(bl_V1locRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for small condition
                                bl_V1loc(bl_V1locRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for small condition

                                bl_V1locRow = bl_V1locRow+1; % advance counter
                            end

                        elseif AK_whichPattern(data.logFilename{iF,iA},contrastStr,true)==1 && istable(data.summaryStats{iF,iA}) % is this a contrast scan?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,contrast(contrastRow,5:61))) % check to advance counter
                                contrastRow = contrastRow+1; % advance counter
                            elseif any(~cellfun(@isempty,contrast(contrastRow,5:61))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(contrastRow) ' of contrast']) % message
                                contrastRow = contrastRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            contrast{contrastRow,1} = subjects{iS}; % store subject name
                            contrast{contrastRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            contrast{contrastRow,3} = data.setN{iF,iA}; % store set number
                            contrast{contrastRow,4} = num2str(scanN); % store scan#

                            contrast(contrastRow,5) = table2cell(data.summaryStats{iF,iA}('fix','noTrackTime')); % store no tracking time for fix condition
                            contrast(contrastRow,6) = table2cell(data.summaryStats{iF,iA}('fix','fixationTime')); % store fixation time for fix condition
                            contrast(contrastRow,7) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedFixationTime')); % store drift-corrected fixation time for fix condition
                            contrast(contrastRow,8) = table2cell(data.summaryStats{iF,iA}('fix','saccadeTime')); % store saccade and blink time for fix condition
                            contrast(contrastRow,9) = table2cell(data.summaryStats{iF,iA}('fix','distanceMean')); % store distance M for fix condition
                            contrast(contrastRow,10) = table2cell(data.summaryStats{iF,iA}('fix','distanceSD')); % store distance SD for fix condition
                            contrast(contrastRow,11) = table2cell(data.summaryStats{iF,iA}('fix','angleMean')); % store angle M for fix condition
                            contrast(contrastRow,12) = table2cell(data.summaryStats{iF,iA}('fix','angleSD')); % store angle SD for fix condition
                            contrast(contrastRow,13) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedDistanceMean')); % store drift-corrected distance M for fix condition
                            contrast(contrastRow,14) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedDistanceSD')); % store drift-corrected distance SD for fix condition
                            contrast(contrastRow,15) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedAngleMean')); % store drift-corrected angle M for fix condition
                            contrast(contrastRow,16) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedAngleSD')); % store drift-corrected angle SD for fix condition
                            contrast(contrastRow,17) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeN')); % store number of saccades for fix condition
                            contrast(contrastRow,18) = table2cell(data.summaryStats{iF,iA}('fix','bigSaccadeN')); % store number of big saccades for fix condition
                            contrast(contrastRow,19) = table2cell(data.summaryStats{iF,iA}('fix','BlinkN')); % store number of saccades for fix condition
                            contrast(contrastRow,20) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeAmplitudeMean')); % store saccade amplitude mean for fix condition
                            contrast(contrastRow,21) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeAmplitudeSD')); % store saccade amplitude SD for fix condition
                            contrast(contrastRow,22) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeVelocityMean')); % store saccade velocity mean for fix condition
                            contrast(contrastRow,23) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeVelocitySD')); % store saccade velocity SD for fix condition

                            contrast(contrastRow,24) = table2cell(data.summaryStats{iF,iA}('low','noTrackTime')); % store no tracking time for low condition
                            contrast(contrastRow,25) = table2cell(data.summaryStats{iF,iA}('low','fixationTime')); % store fixation time for low condition
                            contrast(contrastRow,26) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedFixationTime')); % store drift-corrected fixation time for low condition
                            contrast(contrastRow,27) = table2cell(data.summaryStats{iF,iA}('low','saccadeTime')); % store saccade and blink time for low condition
                            contrast(contrastRow,28) = table2cell(data.summaryStats{iF,iA}('low','distanceMean')); % store distance M for low condition
                            contrast(contrastRow,29) = table2cell(data.summaryStats{iF,iA}('low','distanceSD')); % store distance SD for low condition
                            contrast(contrastRow,30) = table2cell(data.summaryStats{iF,iA}('low','angleMean')); % store angle M for low condition
                            contrast(contrastRow,31) = table2cell(data.summaryStats{iF,iA}('low','angleSD')); % store angle SD for low condition
                            contrast(contrastRow,32) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedDistanceMean')); % store drift-corrected distance M for low condition
                            contrast(contrastRow,33) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedDistanceSD')); % store drift-corrected distance SD for low condition
                            contrast(contrastRow,34) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedAngleMean')); % store drift-corrected angle M for low condition
                            contrast(contrastRow,35) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedAngleSD')); % store drift-corrected angle SD for low condition
                            contrast(contrastRow,36) = table2cell(data.summaryStats{iF,iA}('low','SaccadeN')); % store number of saccades for low condition
                            contrast(contrastRow,37) = table2cell(data.summaryStats{iF,iA}('low','bigSaccadeN')); % store number of big saccades for low condition
                            contrast(contrastRow,38) = table2cell(data.summaryStats{iF,iA}('low','BlinkN')); % store number of saccades for low condition
                            contrast(contrastRow,39) = table2cell(data.summaryStats{iF,iA}('low','SaccadeAmplitudeMean')); % store saccade amplitude mean for low condition
                            contrast(contrastRow,40) = table2cell(data.summaryStats{iF,iA}('low','SaccadeAmplitudeSD')); % store saccade amplitude SD for low condition
                            contrast(contrastRow,41) = table2cell(data.summaryStats{iF,iA}('low','SaccadeVelocityMean')); % store saccade velocity mean for low condition
                            contrast(contrastRow,42) = table2cell(data.summaryStats{iF,iA}('low','SaccadeVelocitySD')); % store saccade velocity SD for low condition

                            contrast(contrastRow,43) = table2cell(data.summaryStats{iF,iA}('hi','noTrackTime')); % store no tracking time for hi condition
                            contrast(contrastRow,44) = table2cell(data.summaryStats{iF,iA}('hi','fixationTime')); % store fixation time for hi condition
                            contrast(contrastRow,45) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedFixationTime')); % store drift-corrected fixation time for hi condition
                            contrast(contrastRow,46) = table2cell(data.summaryStats{iF,iA}('hi','saccadeTime')); % store saccade and blink time for hi condition
                            contrast(contrastRow,47) = table2cell(data.summaryStats{iF,iA}('hi','distanceMean')); % store distance M for hi condition
                            contrast(contrastRow,48) = table2cell(data.summaryStats{iF,iA}('hi','distanceSD')); % store distance SD for hi condition
                            contrast(contrastRow,49) = table2cell(data.summaryStats{iF,iA}('hi','angleMean')); % store angle M for hi condition
                            contrast(contrastRow,50) = table2cell(data.summaryStats{iF,iA}('hi','angleSD')); % store angle SD for hi condition
                            contrast(contrastRow,51) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedDistanceMean')); % store drift-corrected distance M for hi condition
                            contrast(contrastRow,52) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedDistanceSD')); % store drift-corrected distance SD for hi condition
                            contrast(contrastRow,53) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedAngleMean')); % store drift-corrected angle M for hi condition
                            contrast(contrastRow,54) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedAngleSD')); % store drift-corrected angle SD for hi condition
                            contrast(contrastRow,55) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeN')); % store number of saccades for hi condition
                            contrast(contrastRow,56) = table2cell(data.summaryStats{iF,iA}('hi','bigSaccadeN')); % store number of big saccades for hi condition
                            contrast(contrastRow,57) = table2cell(data.summaryStats{iF,iA}('hi','BlinkN')); % store number of saccades for hi condition
                            contrast(contrastRow,58) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeAmplitudeMean')); % store saccade amplitude mean for hi condition
                            contrast(contrastRow,59) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeAmplitudeSD')); % store saccade amplitude SD for hi condition
                            contrast(contrastRow,60) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeVelocityMean')); % store saccade velocity mean for hi condition
                            contrast(contrastRow,61) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeVelocitySD')); % store saccade velocity SD for hi condition

                            contrast(contrastRow,62) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            contrast(contrastRow,63) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            contrast(contrastRow,64) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            contrast(contrastRow,65) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            contrast(contrastRow,66) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            contrast(contrastRow,67) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            contrast(contrastRow,68) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            contrast(contrastRow,69) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            contrast(contrastRow,70) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            contrast(contrastRow,71) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            contrast(contrastRow,72) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            contrast(contrastRow,73) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            contrast(contrastRow,74) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            contrast(contrastRow,75) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            contrast(contrastRow,76) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            contrast(contrastRow,77) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            contrast(contrastRow,78) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            contrast(contrastRow,79) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            contrast(contrastRow,80) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_contrast{bl_contrastRow,1} = subjects{iS}; % store subject name
                                bl_contrast{bl_contrastRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_contrast{bl_contrastRow,3} = data.setN{iF,iA}; % store set number
                                bl_contrast{bl_contrastRow,4} = num2str(scanN); % store scan#
                                bl_contrast(bl_contrastRow,5) = {iB}; % store block#

                                bl_contrast(bl_contrastRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for fix condition
                                bl_contrast(bl_contrastRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for fix condition
                                bl_contrast(bl_contrastRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for fix condition
                                bl_contrast(bl_contrastRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for fix condition
                                bl_contrast(bl_contrastRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''distanceMean''}{1}(iB)}','cell'); % store distance M for fix condition
                                bl_contrast(bl_contrastRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for fix condition
                                bl_contrast(bl_contrastRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''angleMean''}{1}(iB)}','cell'); % store angle M for fix condition
                                bl_contrast(bl_contrastRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''angleSD''}{1}(iB)}','cell'); % store angle SD for fix condition
                                bl_contrast(bl_contrastRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for fix condition
                                bl_contrast(bl_contrastRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for fix condition
                                bl_contrast(bl_contrastRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for fix condition
                                bl_contrast(bl_contrastRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for fix condition
                                bl_contrast(bl_contrastRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for fix condition
                                bl_contrast(bl_contrastRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for fix condition
                                bl_contrast(bl_contrastRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for fix condition
                                bl_contrast(bl_contrastRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for fix condition
                                bl_contrast(bl_contrastRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for fix condition
                                bl_contrast(bl_contrastRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for fix condition
                                bl_contrast(bl_contrastRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for fix condition

                                bl_contrast(bl_contrastRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for low condition
                                bl_contrast(bl_contrastRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for low condition
                                bl_contrast(bl_contrastRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for low condition
                                bl_contrast(bl_contrastRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for low condition
                                bl_contrast(bl_contrastRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''distanceMean''}{1}(iB)}','cell'); % store distance M for low condition
                                bl_contrast(bl_contrastRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for low condition
                                bl_contrast(bl_contrastRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''angleMean''}{1}(iB)}','cell'); % store angle M for low condition
                                bl_contrast(bl_contrastRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''angleSD''}{1}(iB)}','cell'); % store angle SD for low condition
                                bl_contrast(bl_contrastRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for low condition
                                bl_contrast(bl_contrastRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for low condition
                                bl_contrast(bl_contrastRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for low condition
                                bl_contrast(bl_contrastRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for low condition
                                bl_contrast(bl_contrastRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for low condition
                                bl_contrast(bl_contrastRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for low condition
                                bl_contrast(bl_contrastRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for low condition
                                bl_contrast(bl_contrastRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for low condition
                                bl_contrast(bl_contrastRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for low condition
                                bl_contrast(bl_contrastRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for low condition
                                bl_contrast(bl_contrastRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for low condition

                                bl_contrast(bl_contrastRow,44) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for hi condition
                                bl_contrast(bl_contrastRow,45) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for hi condition
                                bl_contrast(bl_contrastRow,46) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for hi condition
                                bl_contrast(bl_contrastRow,47) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for hi condition
                                bl_contrast(bl_contrastRow,48) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''distanceMean''}{1}(iB)}','cell'); % store distance M for hi condition
                                bl_contrast(bl_contrastRow,49) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for hi condition
                                bl_contrast(bl_contrastRow,50) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''angleMean''}{1}(iB)}','cell'); % store angle M for hi condition
                                bl_contrast(bl_contrastRow,51) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''angleSD''}{1}(iB)}','cell'); % store angle SD for hi condition
                                bl_contrast(bl_contrastRow,52) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for hi condition
                                bl_contrast(bl_contrastRow,53) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for hi condition
                                bl_contrast(bl_contrastRow,54) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for hi condition
                                bl_contrast(bl_contrastRow,55) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for hi condition
                                bl_contrast(bl_contrastRow,56) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for hi condition
                                bl_contrast(bl_contrastRow,57) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for hi condition
                                bl_contrast(bl_contrastRow,58) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for hi condition
                                bl_contrast(bl_contrastRow,59) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for hi condition
                                bl_contrast(bl_contrastRow,60) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for hi condition
                                bl_contrast(bl_contrastRow,61) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for hi condition
                                bl_contrast(bl_contrastRow,62) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for hi condition

                                bl_contrastRow = bl_contrastRow+1; % advance counter
                            end
                            
                        elseif AK_whichPattern(data.logFilename{iF,iA},supStr,true)==1 && istable(data.summaryStats{iF,iA}) % is this a suppression scan?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,sup(supRow,5:42))) % check to advance counter
                                supRow = supRow+1; % advance counter
                            elseif any(~cellfun(@isempty,sup(supRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(supRow) ' of sup']) % message
                                supRow = supRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            sup{supRow,1} = subjects{iS}; % store subject name
                            sup{supRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            sup{supRow,3} = data.setN{iF,iA}; % store set number
                            sup{supRow,4} = num2str(scanN); % store scan#

                            sup(supRow,5) = table2cell(data.summaryStats{iF,iA}('small','noTrackTime')); % store no tracking time for small condition
                            sup(supRow,6) = table2cell(data.summaryStats{iF,iA}('small','fixationTime')); % store fixation time for small condition
                            sup(supRow,7) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedFixationTime')); % store drift-corrected fixation time for small condition
                            sup(supRow,8) = table2cell(data.summaryStats{iF,iA}('small','saccadeTime')); % store saccade and blink time for small condition
                            sup(supRow,9) = table2cell(data.summaryStats{iF,iA}('small','distanceMean')); % store distance M for small condition
                            sup(supRow,10) = table2cell(data.summaryStats{iF,iA}('small','distanceSD')); % store distance SD for small condition
                            sup(supRow,11) = table2cell(data.summaryStats{iF,iA}('small','angleMean')); % store angle M for small condition
                            sup(supRow,12) = table2cell(data.summaryStats{iF,iA}('small','angleSD')); % store angle SD for small condition
                            sup(supRow,13) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedDistanceMean')); % store drift-corrected distance M for small condition
                            sup(supRow,14) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedDistanceSD')); % store drift-corrected distance SD for small condition
                            sup(supRow,15) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedAngleMean')); % store drift-corrected angle M for small condition
                            sup(supRow,16) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedAngleSD')); % store drift-corrected angle SD for small condition
                            sup(supRow,17) = table2cell(data.summaryStats{iF,iA}('small','SaccadeN')); % store number of saccades for small condition
                            sup(supRow,18) = table2cell(data.summaryStats{iF,iA}('small','bigSaccadeN')); % store number of big saccades for small condition
                            sup(supRow,19) = table2cell(data.summaryStats{iF,iA}('small','BlinkN')); % store number of saccades for small condition
                            sup(supRow,20) = table2cell(data.summaryStats{iF,iA}('small','SaccadeAmplitudeMean')); % store saccade amplitude mean for small condition
                            sup(supRow,21) = table2cell(data.summaryStats{iF,iA}('small','SaccadeAmplitudeSD')); % store saccade amplitude SD for small condition
                            sup(supRow,22) = table2cell(data.summaryStats{iF,iA}('small','SaccadeVelocityMean')); % store saccade velocity mean for small condition
                            sup(supRow,23) = table2cell(data.summaryStats{iF,iA}('small','SaccadeVelocitySD')); % store saccade velocity SD for small condition

                            sup(supRow,24) = table2cell(data.summaryStats{iF,iA}('big','noTrackTime')); % store no tracking time for big condition
                            sup(supRow,25) = table2cell(data.summaryStats{iF,iA}('big','fixationTime')); % store fixation time for big condition
                            sup(supRow,26) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedFixationTime')); % store drift-corrected fixation time for big condition
                            sup(supRow,27) = table2cell(data.summaryStats{iF,iA}('big','saccadeTime')); % store saccade and blink time for big condition
                            sup(supRow,28) = table2cell(data.summaryStats{iF,iA}('big','distanceMean')); % store distance M for big condition
                            sup(supRow,29) = table2cell(data.summaryStats{iF,iA}('big','distanceSD')); % store distance SD for big condition
                            sup(supRow,30) = table2cell(data.summaryStats{iF,iA}('big','angleMean')); % store angle M for big condition
                            sup(supRow,31) = table2cell(data.summaryStats{iF,iA}('big','angleSD')); % store angle SD for big condition
                            sup(supRow,32) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedDistanceMean')); % store drift-corrected distance M for big condition
                            sup(supRow,33) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedDistanceSD')); % store drift-corrected distance SD for big condition
                            sup(supRow,34) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedAngleMean')); % store drift-corrected angle M for big condition
                            sup(supRow,35) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedAngleSD')); % store drift-corrected angle SD for big condition
                            sup(supRow,36) = table2cell(data.summaryStats{iF,iA}('big','SaccadeN')); % store number of saccades for big condition
                            sup(supRow,37) = table2cell(data.summaryStats{iF,iA}('big','bigSaccadeN')); % store number of big saccades for big condition
                            sup(supRow,38) = table2cell(data.summaryStats{iF,iA}('big','BlinkN')); % store number of saccades for big condition
                            sup(supRow,39) = table2cell(data.summaryStats{iF,iA}('big','SaccadeAmplitudeMean')); % store saccade amplitude mean for big condition
                            sup(supRow,40) = table2cell(data.summaryStats{iF,iA}('big','SaccadeAmplitudeSD')); % store saccade amplitude SD for big condition
                            sup(supRow,41) = table2cell(data.summaryStats{iF,iA}('big','SaccadeVelocityMean')); % store saccade velocity mean for big condition
                            sup(supRow,42) = table2cell(data.summaryStats{iF,iA}('big','SaccadeVelocitySD')); % store saccade velocity SD for big condition

                            sup(supRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            sup(supRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            sup(supRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            sup(supRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            sup(supRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            sup(supRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            sup(supRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            sup(supRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            sup(supRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            sup(supRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            sup(supRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            sup(supRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            sup(supRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            sup(supRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            sup(supRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            sup(supRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            sup(supRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            sup(supRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            sup(supRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_sup{bl_supRow,1} = subjects{iS}; % store subject name
                                bl_sup{bl_supRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_sup{bl_supRow,3} = data.setN{iF,iA}; % store set number 
                                bl_sup{bl_supRow,4} = num2str(scanN); % store scan#
                                bl_sup(bl_supRow,5) = {iB}; % store block#

                                bl_sup(bl_supRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for small condition
                                bl_sup(bl_supRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for small condition
                                bl_sup(bl_supRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for small condition
                                bl_sup(bl_supRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for small condition
                                bl_sup(bl_supRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''distanceMean''}{1}(iB)}','cell'); % store distance M for small condition
                                bl_sup(bl_supRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for small condition
                                bl_sup(bl_supRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''angleMean''}{1}(iB)}','cell'); % store angle M for small condition
                                bl_sup(bl_supRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''angleSD''}{1}(iB)}','cell'); % store angle SD for small condition
                                bl_sup(bl_supRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for small condition
                                bl_sup(bl_supRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for small condition
                                bl_sup(bl_supRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for small condition
                                bl_sup(bl_supRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for small condition
                                bl_sup(bl_supRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for small condition
                                bl_sup(bl_supRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for small condition
                                bl_sup(bl_supRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for small condition
                                bl_sup(bl_supRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for small condition
                                bl_sup(bl_supRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for small condition
                                bl_sup(bl_supRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for small condition
                                bl_sup(bl_supRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for small condition

                                bl_sup(bl_supRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for big condition
                                bl_sup(bl_supRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for big condition
                                bl_sup(bl_supRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for big condition
                                bl_sup(bl_supRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for big condition
                                bl_sup(bl_supRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''distanceMean''}{1}(iB)}','cell'); % store distance M for big condition
                                bl_sup(bl_supRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for big condition
                                bl_sup(bl_supRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''angleMean''}{1}(iB)}','cell'); % store angle M for big condition
                                bl_sup(bl_supRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''angleSD''}{1}(iB)}','cell'); % store angle SD for big condition
                                bl_sup(bl_supRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for big condition
                                bl_sup(bl_supRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for big condition
                                bl_sup(bl_supRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for big condition
                                bl_sup(bl_supRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for big condition
                                bl_sup(bl_supRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for big condition
                                bl_sup(bl_supRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for big condition
                                bl_sup(bl_supRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for big condition
                                bl_sup(bl_supRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for big condition
                                bl_sup(bl_supRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for big condition
                                bl_sup(bl_supRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for big condition
                                bl_sup(bl_supRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for big condition

                                bl_supRow = bl_supRow+1; % advance counter
                            end

                        elseif AK_whichPattern(data.logFilename{iF,iA},sumStr,true)==1 && istable(data.summaryStats{iF,iA}) % is this a summation scan?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,sum(sumRow,5:42))) % check to advance counter
                                sumRow = sumRow+1; % advance counter
                            elseif any(~cellfun(@isempty,sum(sumRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(sumRow) ' of sum']) % message
                                sumRow = sumRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            sum{sumRow,1} = subjects{iS}; % store subject name
                            sum{sumRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            sum{sumRow,3} = data.setN{iF,iA}; % store set number
                            sum{sumRow,4} = num2str(scanN); % store scan#

                            sum(sumRow,5) = table2cell(data.summaryStats{iF,iA}('small','noTrackTime')); % store no tracking time for small condition
                            sum(sumRow,6) = table2cell(data.summaryStats{iF,iA}('small','fixationTime')); % store fixation time for small condition
                            sum(sumRow,7) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedFixationTime')); % store drift-corrected fixation time for small condition
                            sum(sumRow,8) = table2cell(data.summaryStats{iF,iA}('small','saccadeTime')); % store saccade and blink time for small condition
                            sum(sumRow,9) = table2cell(data.summaryStats{iF,iA}('small','distanceMean')); % store distance M for small condition
                            sum(sumRow,10) = table2cell(data.summaryStats{iF,iA}('small','distanceSD')); % store distance SD for small condition
                            sum(sumRow,11) = table2cell(data.summaryStats{iF,iA}('small','angleMean')); % store angle M for small condition
                            sum(sumRow,12) = table2cell(data.summaryStats{iF,iA}('small','angleSD')); % store angle SD for small condition
                            sum(sumRow,13) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedDistanceMean')); % store drift-corrected distance M for small condition
                            sum(sumRow,14) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedDistanceSD')); % store drift-corrected distance SD for small condition
                            sum(sumRow,15) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedAngleMean')); % store drift-corrected angle M for small condition
                            sum(sumRow,16) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedAngleSD')); % store drift-corrected angle SD for small condition
                            sum(sumRow,17) = table2cell(data.summaryStats{iF,iA}('small','SaccadeN')); % store number of saccades for small condition
                            sum(sumRow,18) = table2cell(data.summaryStats{iF,iA}('small','bigSaccadeN')); % store number of big saccades for small condition
                            sum(sumRow,19) = table2cell(data.summaryStats{iF,iA}('small','BlinkN')); % store number of saccades for small condition
                            sum(sumRow,20) = table2cell(data.summaryStats{iF,iA}('small','SaccadeAmplitudeMean')); % store saccade amplitude mean for small condition
                            sum(sumRow,21) = table2cell(data.summaryStats{iF,iA}('small','SaccadeAmplitudeSD')); % store saccade amplitude SD for small condition
                            sum(sumRow,22) = table2cell(data.summaryStats{iF,iA}('small','SaccadeVelocityMean')); % store saccade velocity mean for small condition
                            sum(sumRow,23) = table2cell(data.summaryStats{iF,iA}('small','SaccadeVelocitySD')); % store saccade velocity SD for small condition

                            sum(sumRow,24) = table2cell(data.summaryStats{iF,iA}('big','noTrackTime')); % store no tracking time for big condition
                            sum(sumRow,25) = table2cell(data.summaryStats{iF,iA}('big','fixationTime')); % store fixation time for big condition
                            sum(sumRow,26) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedFixationTime')); % store drift-corrected fixation time for big condition
                            sum(sumRow,27) = table2cell(data.summaryStats{iF,iA}('big','saccadeTime')); % store saccade and blink time for big condition
                            sum(sumRow,28) = table2cell(data.summaryStats{iF,iA}('big','distanceMean')); % store distance M for big condition
                            sum(sumRow,29) = table2cell(data.summaryStats{iF,iA}('big','distanceSD')); % store distance SD for big condition
                            sum(sumRow,30) = table2cell(data.summaryStats{iF,iA}('big','angleMean')); % store angle M for big condition
                            sum(sumRow,31) = table2cell(data.summaryStats{iF,iA}('big','angleSD')); % store angle SD for big condition
                            sum(sumRow,32) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedDistanceMean')); % store drift-corrected distance M for big condition
                            sum(sumRow,33) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedDistanceSD')); % store drift-corrected distance SD for big condition
                            sum(sumRow,34) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedAngleMean')); % store drift-corrected angle M for big condition
                            sum(sumRow,35) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedAngleSD')); % store drift-corrected angle SD for big condition
                            sum(sumRow,36) = table2cell(data.summaryStats{iF,iA}('big','SaccadeN')); % store number of saccades for big condition
                            sum(sumRow,37) = table2cell(data.summaryStats{iF,iA}('big','bigSaccadeN')); % store number of big saccades for big condition
                            sum(sumRow,38) = table2cell(data.summaryStats{iF,iA}('big','BlinkN')); % store number of saccades for big condition
                            sum(sumRow,39) = table2cell(data.summaryStats{iF,iA}('big','SaccadeAmplitudeMean')); % store saccade amplitude mean for big condition
                            sum(sumRow,40) = table2cell(data.summaryStats{iF,iA}('big','SaccadeAmplitudeSD')); % store saccade amplitude SD for big condition
                            sum(sumRow,41) = table2cell(data.summaryStats{iF,iA}('big','SaccadeVelocityMean')); % store saccade velocity mean for big condition
                            sum(sumRow,42) = table2cell(data.summaryStats{iF,iA}('big','SaccadeVelocitySD')); % store saccade velocity SD for big condition

                            sum(sumRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            sum(sumRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            sum(sumRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            sum(sumRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            sum(sumRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            sum(sumRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            sum(sumRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            sum(sumRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            sum(sumRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            sum(sumRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            sum(sumRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            sum(sumRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            sum(sumRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            sum(sumRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            sum(sumRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            sum(sumRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            sum(sumRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            sum(sumRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            sum(sumRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_sum{bl_sumRow,1} = subjects{iS}; % store subject name
                                bl_sum{bl_sumRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_sum{bl_sumRow,3} = data.setN{iF,iA}; % store set number
                                bl_sum{bl_sumRow,4} = num2str(scanN); % store scan#
                                bl_sum(bl_sumRow,5) = {iB}; % store block#

                                bl_sum(bl_sumRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for small condition
                                bl_sum(bl_sumRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for small condition
                                bl_sum(bl_sumRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for small condition
                                bl_sum(bl_sumRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for small condition
                                bl_sum(bl_sumRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''distanceMean''}{1}(iB)}','cell'); % store distance M for small condition
                                bl_sum(bl_sumRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for small condition
                                bl_sum(bl_sumRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''angleMean''}{1}(iB)}','cell'); % store angle M for small condition
                                bl_sum(bl_sumRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''angleSD''}{1}(iB)}','cell'); % store angle SD for small condition
                                bl_sum(bl_sumRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for small condition
                                bl_sum(bl_sumRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for small condition
                                bl_sum(bl_sumRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for small condition
                                bl_sum(bl_sumRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for small condition
                                bl_sum(bl_sumRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for small condition
                                bl_sum(bl_sumRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for small condition
                                bl_sum(bl_sumRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for small condition
                                bl_sum(bl_sumRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for small condition
                                bl_sum(bl_sumRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for small condition
                                bl_sum(bl_sumRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for small condition
                                bl_sum(bl_sumRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for small condition

                                bl_sum(bl_sumRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for big condition
                                bl_sum(bl_sumRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for big condition
                                bl_sum(bl_sumRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for big condition
                                bl_sum(bl_sumRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for big condition
                                bl_sum(bl_sumRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''distanceMean''}{1}(iB)}','cell'); % store distance M for big condition
                                bl_sum(bl_sumRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for big condition
                                bl_sum(bl_sumRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''angleMean''}{1}(iB)}','cell'); % store angle M for big condition
                                bl_sum(bl_sumRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''angleSD''}{1}(iB)}','cell'); % store angle SD for big condition
                                bl_sum(bl_sumRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for big condition
                                bl_sum(bl_sumRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for big condition
                                bl_sum(bl_sumRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for big condition
                                bl_sum(bl_sumRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for big condition
                                bl_sum(bl_sumRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for big condition
                                bl_sum(bl_sumRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for big condition
                                bl_sum(bl_sumRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for big condition
                                bl_sum(bl_sumRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for big condition
                                bl_sum(bl_sumRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for big condition
                                bl_sum(bl_sumRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for big condition
                                bl_sum(bl_sumRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for big condition

                                bl_sumRow = bl_sumRow+1; % advance counter
                            end
                            
                        elseif AK_whichPattern(data.logFilename{iF,iA},ftapStr,true)==1 && istable(data.summaryStats{iF,iA}) % is this a ftap scan?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,ftap(ftapRow,5:61))) % check to advance counter
                                ftapRow = ftapRow+1; % advance counter
                            elseif any(~cellfun(@isempty,ftap(ftapRow,5:61))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(ftapRow) ' of ftap']) % message
                                ftapRow = ftapRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            ftap{ftapRow,1} = subjects{iS}; % store subject name
                            ftap{ftapRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            ftap{ftapRow,3} = data.setN{iF,iA}; % store set number
                            ftap{ftapRow,4} = num2str(scanN); % store scan#

                            ftap(ftapRow,5) = table2cell(data.summaryStats{iF,iA}('fix','noTrackTime')); % store no tracking time for fix condition
                            ftap(ftapRow,6) = table2cell(data.summaryStats{iF,iA}('fix','fixationTime')); % store fixation time for fix condition
                            ftap(ftapRow,7) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedFixationTime')); % store drift-corrected fixation time for fix condition
                            ftap(ftapRow,8) = table2cell(data.summaryStats{iF,iA}('fix','saccadeTime')); % store saccade and blink time for fix condition
                            ftap(ftapRow,9) = table2cell(data.summaryStats{iF,iA}('fix','distanceMean')); % store distance M for fix condition
                            ftap(ftapRow,10) = table2cell(data.summaryStats{iF,iA}('fix','distanceSD')); % store distance SD for fix condition
                            ftap(ftapRow,11) = table2cell(data.summaryStats{iF,iA}('fix','angleMean')); % store angle M for fix condition
                            ftap(ftapRow,12) = table2cell(data.summaryStats{iF,iA}('fix','angleSD')); % store angle SD for fix condition
                            ftap(ftapRow,13) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedDistanceMean')); % store drift-corrected distance M for fix condition
                            ftap(ftapRow,14) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedDistanceSD')); % store drift-corrected distance SD for fix condition
                            ftap(ftapRow,15) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedAngleMean')); % store drift-corrected angle M for fix condition
                            ftap(ftapRow,16) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedAngleSD')); % store drift-corrected angle SD for fix condition
                            ftap(ftapRow,17) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeN')); % store number of saccades for fix condition
                            ftap(ftapRow,18) = table2cell(data.summaryStats{iF,iA}('fix','bigSaccadeN')); % store number of big saccades for fix condition
                            ftap(ftapRow,19) = table2cell(data.summaryStats{iF,iA}('fix','BlinkN')); % store number of saccades for fix condition
                            ftap(ftapRow,20) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeAmplitudeMean')); % store saccade amplitude mean for fix condition
                            ftap(ftapRow,21) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeAmplitudeSD')); % store saccade amplitude SD for fix condition
                            ftap(ftapRow,22) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeVelocityMean')); % store saccade velocity mean for fix condition
                            ftap(ftapRow,23) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeVelocitySD')); % store saccade velocity SD for fix condition

                            ftap(ftapRow,24) = table2cell(data.summaryStats{iF,iA}('standard','noTrackTime')); % store no tracking time for predictable condition
                            ftap(ftapRow,25) = table2cell(data.summaryStats{iF,iA}('standard','fixationTime')); % store fixation time for predictable condition
                            ftap(ftapRow,26) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedFixationTime')); % store drift-corrected fixation time for predictable condition
                            ftap(ftapRow,27) = table2cell(data.summaryStats{iF,iA}('standard','saccadeTime')); % store saccade and blink time for predictable condition
                            ftap(ftapRow,28) = table2cell(data.summaryStats{iF,iA}('standard','distanceMean')); % store distance M for predictable condition
                            ftap(ftapRow,29) = table2cell(data.summaryStats{iF,iA}('standard','distanceSD')); % store distance SD for predictable condition
                            ftap(ftapRow,30) = table2cell(data.summaryStats{iF,iA}('standard','angleMean')); % store angle M for predictable condition
                            ftap(ftapRow,31) = table2cell(data.summaryStats{iF,iA}('standard','angleSD')); % store angle SD for predictable condition
                            ftap(ftapRow,32) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedDistanceMean')); % store drift-corrected distance M for predictable condition
                            ftap(ftapRow,33) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedDistanceSD')); % store drift-corrected distance SD for predictable condition
                            ftap(ftapRow,34) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedAngleMean')); % store drift-corrected angle M for predictable condition
                            ftap(ftapRow,35) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedAngleSD')); % store drift-corrected angle SD for predictable condition
                            ftap(ftapRow,36) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeN')); % store number of saccades for predictable condition
                            ftap(ftapRow,37) = table2cell(data.summaryStats{iF,iA}('standard','bigSaccadeN')); % store number of big saccades for predictable condition
                            ftap(ftapRow,38) = table2cell(data.summaryStats{iF,iA}('standard','BlinkN')); % store number of saccades for predictable condition
                            ftap(ftapRow,39) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeAmplitudeMean')); % store saccade amplitude mean for predictable condition
                            ftap(ftapRow,40) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeAmplitudeSD')); % store saccade amplitude SD for predictable condition
                            ftap(ftapRow,41) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeVelocityMean')); % store saccade velocity mean for predictable condition
                            ftap(ftapRow,42) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeVelocitySD')); % store saccade velocity SD for predictable condition

                            ftap(ftapRow,43) = table2cell(data.summaryStats{iF,iA}('variable','noTrackTime')); % store no tracking time for variable condition
                            ftap(ftapRow,44) = table2cell(data.summaryStats{iF,iA}('variable','fixationTime')); % store fixation time for variable condition
                            ftap(ftapRow,45) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedFixationTime')); % store drift-corrected fixation time for variable condition
                            ftap(ftapRow,46) = table2cell(data.summaryStats{iF,iA}('variable','saccadeTime')); % store saccade and blink time for variable condition
                            ftap(ftapRow,47) = table2cell(data.summaryStats{iF,iA}('variable','distanceMean')); % store distance M for variable condition
                            ftap(ftapRow,48) = table2cell(data.summaryStats{iF,iA}('variable','distanceSD')); % store distance SD for variable condition
                            ftap(ftapRow,49) = table2cell(data.summaryStats{iF,iA}('variable','angleMean')); % store angle M for variable condition
                            ftap(ftapRow,50) = table2cell(data.summaryStats{iF,iA}('variable','angleSD')); % store angle SD for variable condition
                            ftap(ftapRow,51) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedDistanceMean')); % store drift-corrected distance M for variable condition
                            ftap(ftapRow,52) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedDistanceSD')); % store drift-corrected distance SD for variable condition
                            ftap(ftapRow,53) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedAngleMean')); % store drift-corrected angle M for variable condition
                            ftap(ftapRow,54) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedAngleSD')); % store drift-corrected angle SD for variable condition
                            ftap(ftapRow,55) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeN')); % store number of saccades for variable condition
                            ftap(ftapRow,56) = table2cell(data.summaryStats{iF,iA}('variable','bigSaccadeN')); % store number of big saccades for variable condition
                            ftap(ftapRow,57) = table2cell(data.summaryStats{iF,iA}('variable','BlinkN')); % store number of saccades for variable condition
                            ftap(ftapRow,58) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeAmplitudeMean')); % store saccade amplitude mean for variable condition
                            ftap(ftapRow,59) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeAmplitudeSD')); % store saccade amplitude SD for variable condition
                            ftap(ftapRow,60) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeVelocityMean')); % store saccade velocity mean for variable condition
                            ftap(ftapRow,61) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeVelocitySD')); % store saccade velocity SD for variable condition

                            ftap(ftapRow,62) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            ftap(ftapRow,63) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            ftap(ftapRow,64) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            ftap(ftapRow,65) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            ftap(ftapRow,66) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            ftap(ftapRow,67) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            ftap(ftapRow,68) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            ftap(ftapRow,69) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            ftap(ftapRow,70) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            ftap(ftapRow,71) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            ftap(ftapRow,72) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            ftap(ftapRow,73) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            ftap(ftapRow,74) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            ftap(ftapRow,75) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            ftap(ftapRow,76) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            ftap(ftapRow,77) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            ftap(ftapRow,78) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            ftap(ftapRow,79) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            ftap(ftapRow,80) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_ftap{bl_ftapRow,1} = subjects{iS}; % store subject name
                                bl_ftap{bl_ftapRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_ftap{bl_ftapRow,3} = data.setN{iF,iA}; % store set number
                                bl_ftap{bl_ftapRow,4} = num2str(scanN); % store scan#
                                bl_ftap(bl_ftapRow,5) = {iB}; % store block#

                                bl_ftap(bl_ftapRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for fix condition
                                bl_ftap(bl_ftapRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for fix condition
                                bl_ftap(bl_ftapRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for fix condition
                                bl_ftap(bl_ftapRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for fix condition
                                bl_ftap(bl_ftapRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''distanceMean''}{1}(iB)}','cell'); % store distance M for fix condition
                                bl_ftap(bl_ftapRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for fix condition
                                bl_ftap(bl_ftapRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''angleMean''}{1}(iB)}','cell'); % store angle M for fix condition
                                bl_ftap(bl_ftapRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''angleSD''}{1}(iB)}','cell'); % store angle SD for fix condition
                                bl_ftap(bl_ftapRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for fix condition
                                bl_ftap(bl_ftapRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for fix condition
                                bl_ftap(bl_ftapRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for fix condition
                                bl_ftap(bl_ftapRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for fix condition
                                bl_ftap(bl_ftapRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for fix condition
                                bl_ftap(bl_ftapRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for fix condition
                                bl_ftap(bl_ftapRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for fix condition
                                bl_ftap(bl_ftapRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for fix condition
                                bl_ftap(bl_ftapRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for fix condition
                                bl_ftap(bl_ftapRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for fix condition
                                bl_ftap(bl_ftapRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for fix condition

                                bl_ftap(bl_ftapRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for predictable condition
                                bl_ftap(bl_ftapRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for predictable condition
                                bl_ftap(bl_ftapRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for predictable condition
                                bl_ftap(bl_ftapRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for predictable condition
                                bl_ftap(bl_ftapRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''distanceMean''}{1}(iB)}','cell'); % store distance M for predictable condition
                                bl_ftap(bl_ftapRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for predictable condition
                                bl_ftap(bl_ftapRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''angleMean''}{1}(iB)}','cell'); % store angle M for predictable condition
                                bl_ftap(bl_ftapRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''angleSD''}{1}(iB)}','cell'); % store angle SD for predictable condition
                                bl_ftap(bl_ftapRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for predictable condition
                                bl_ftap(bl_ftapRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for predictable condition
                                bl_ftap(bl_ftapRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for predictable condition
                                bl_ftap(bl_ftapRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for predictable condition
                                bl_ftap(bl_ftapRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for predictable condition
                                bl_ftap(bl_ftapRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for predictable condition
                                bl_ftap(bl_ftapRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for predictable condition
                                bl_ftap(bl_ftapRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for predictable condition
                                bl_ftap(bl_ftapRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for predictable condition
                                bl_ftap(bl_ftapRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for predictable condition
                                bl_ftap(bl_ftapRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for predictable condition

                                bl_ftap(bl_ftapRow,44) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for variable condition
                                bl_ftap(bl_ftapRow,45) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for variable condition
                                bl_ftap(bl_ftapRow,46) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for variable condition
                                bl_ftap(bl_ftapRow,47) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for variable condition
                                bl_ftap(bl_ftapRow,48) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''distanceMean''}{1}(iB)}','cell'); % store distance M for variable condition
                                bl_ftap(bl_ftapRow,49) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for variable condition
                                bl_ftap(bl_ftapRow,50) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''angleMean''}{1}(iB)}','cell'); % store angle M for variable condition
                                bl_ftap(bl_ftapRow,51) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''angleSD''}{1}(iB)}','cell'); % store angle SD for variable condition
                                bl_ftap(bl_ftapRow,52) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for variable condition
                                bl_ftap(bl_ftapRow,53) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for variable condition
                                bl_ftap(bl_ftapRow,54) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for variable condition
                                bl_ftap(bl_ftapRow,55) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for variable condition
                                bl_ftap(bl_ftapRow,56) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for variable condition
                                bl_ftap(bl_ftapRow,57) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for variable condition
                                bl_ftap(bl_ftapRow,58) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for variable condition
                                bl_ftap(bl_ftapRow,59) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for variable condition
                                bl_ftap(bl_ftapRow,60) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for variable condition
                                bl_ftap(bl_ftapRow,61) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for variable condition
                                bl_ftap(bl_ftapRow,62) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for variable condition

                                bl_ftapRow = bl_ftapRow+1; % advance counter
                            end
                        end
                    case 1 % only good quality data
                        % determine which type of scan iA corresponds to
                        if AK_whichPattern(data.logFilename{iF,iA},MTlocStr,true)==1 && istable(data.summaryStats{iF,iA}) && data.quality{iF,iA}==1 % is this an MT localizer?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,MTloc(MTlocRow,5:42))) % check to advance counter
                                MTlocRow = MTlocRow+1; % advance counter
                            elseif any(~cellfun(@isempty,MTloc(MTlocRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(MTlocRow) ' of MTloc']) % message
                                MTlocRow = MTlocRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            MTloc{MTlocRow,1} = subjects{iS}; % store subject name
                            MTloc{MTlocRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            MTloc{MTlocRow,3} = data.setN{iF,iA}; % store set number
                            MTloc{MTlocRow,4} = num2str(scanN); % store scan#

                            MTloc(MTlocRow,5) = table2cell(data.summaryStats{iF,iA}('moving','noTrackTime')); % store no tracking time for moving condition
                            MTloc(MTlocRow,6) = table2cell(data.summaryStats{iF,iA}('moving','fixationTime')); % store fixation time for moving condition
                            MTloc(MTlocRow,7) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedFixationTime')); % store drift-corrected fixation time for moving condition
                            MTloc(MTlocRow,8) = table2cell(data.summaryStats{iF,iA}('moving','saccadeTime')); % store saccade and blink time for moving condition
                            MTloc(MTlocRow,9) = table2cell(data.summaryStats{iF,iA}('moving','distanceMean')); % store distance M for moving condition
                            MTloc(MTlocRow,10) = table2cell(data.summaryStats{iF,iA}('moving','distanceSD')); % store distance SD for moving condition
                            MTloc(MTlocRow,11) = table2cell(data.summaryStats{iF,iA}('moving','angleMean')); % store angle M for moving condition
                            MTloc(MTlocRow,12) = table2cell(data.summaryStats{iF,iA}('moving','angleSD')); % store angle SD for moving condition
                            MTloc(MTlocRow,13) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedDistanceMean')); % store drift-corrected distance M for moving condition
                            MTloc(MTlocRow,14) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedDistanceSD')); % store drift-corrected distance SD for moving condition
                            MTloc(MTlocRow,15) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedAngleMean')); % store drift-corrected angle M for moving condition
                            MTloc(MTlocRow,16) = table2cell(data.summaryStats{iF,iA}('moving','driftCorrectedAngleSD')); % store drift-corrected angle SD for moving condition
                            MTloc(MTlocRow,17) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeN')); % store number of saccades for moving condition
                            MTloc(MTlocRow,18) = table2cell(data.summaryStats{iF,iA}('moving','bigSaccadeN')); % store number of big saccades for moving condition
                            MTloc(MTlocRow,19) = table2cell(data.summaryStats{iF,iA}('moving','BlinkN')); % store number of saccades for moving condition
                            MTloc(MTlocRow,20) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeAmplitudeMean')); % store saccade amplitude mean for moving condition
                            MTloc(MTlocRow,21) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeAmplitudeSD')); % store saccade amplitude SD for moving condition
                            MTloc(MTlocRow,22) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeVelocityMean')); % store saccade velocity mean for moving condition
                            MTloc(MTlocRow,23) = table2cell(data.summaryStats{iF,iA}('moving','SaccadeVelocitySD')); % store saccade velocity SD for moving condition

                            MTloc(MTlocRow,24) = table2cell(data.summaryStats{iF,iA}('static','noTrackTime')); % store no tracking time for static condition
                            MTloc(MTlocRow,25) = table2cell(data.summaryStats{iF,iA}('static','fixationTime')); % store fixation time for static condition
                            MTloc(MTlocRow,26) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedFixationTime')); % store drift-corrected fixation time for static condition
                            MTloc(MTlocRow,27) = table2cell(data.summaryStats{iF,iA}('static','saccadeTime')); % store saccade and blink time for static condition
                            MTloc(MTlocRow,28) = table2cell(data.summaryStats{iF,iA}('static','distanceMean')); % store distance M for static condition
                            MTloc(MTlocRow,29) = table2cell(data.summaryStats{iF,iA}('static','distanceSD')); % store distance SD for static condition
                            MTloc(MTlocRow,30) = table2cell(data.summaryStats{iF,iA}('static','angleMean')); % store angle M for static condition
                            MTloc(MTlocRow,31) = table2cell(data.summaryStats{iF,iA}('static','angleSD')); % store angle SD for static condition
                            MTloc(MTlocRow,32) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedDistanceMean')); % store drift-corrected distance M for static condition
                            MTloc(MTlocRow,33) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedDistanceSD')); % store drift-corrected distance SD for static condition
                            MTloc(MTlocRow,34) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedAngleMean')); % store drift-corrected angle M for static condition
                            MTloc(MTlocRow,35) = table2cell(data.summaryStats{iF,iA}('static','driftCorrectedAngleSD')); % store drift-corrected angle SD for static condition
                            MTloc(MTlocRow,36) = table2cell(data.summaryStats{iF,iA}('static','SaccadeN')); % store number of saccades for static condition
                            MTloc(MTlocRow,37) = table2cell(data.summaryStats{iF,iA}('static','bigSaccadeN')); % store number of big saccades for static condition
                            MTloc(MTlocRow,38) = table2cell(data.summaryStats{iF,iA}('static','BlinkN')); % store number of saccades for static condition
                            MTloc(MTlocRow,39) = table2cell(data.summaryStats{iF,iA}('static','SaccadeAmplitudeMean')); % store saccade amplitude mean for static condition
                            MTloc(MTlocRow,40) = table2cell(data.summaryStats{iF,iA}('static','SaccadeAmplitudeSD')); % store saccade amplitude SD for static condition
                            MTloc(MTlocRow,41) = table2cell(data.summaryStats{iF,iA}('static','SaccadeVelocityMean')); % store saccade velocity mean for static condition
                            MTloc(MTlocRow,42) = table2cell(data.summaryStats{iF,iA}('static','SaccadeVelocitySD')); % store saccade velocity SD for static condition

                            MTloc(MTlocRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            MTloc(MTlocRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            MTloc(MTlocRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            MTloc(MTlocRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            MTloc(MTlocRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            MTloc(MTlocRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            MTloc(MTlocRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            MTloc(MTlocRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            MTloc(MTlocRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            MTloc(MTlocRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            MTloc(MTlocRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            MTloc(MTlocRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            MTloc(MTlocRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            MTloc(MTlocRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            MTloc(MTlocRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            MTloc(MTlocRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            MTloc(MTlocRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            MTloc(MTlocRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            MTloc(MTlocRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_MTloc{bl_MTlocRow,1} = subjects{iS}; % store subject name
                                bl_MTloc{bl_MTlocRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_MTloc{bl_MTlocRow,3} = data.setN{iF,iA}; % store set number
                                bl_MTloc{bl_MTlocRow,4} = num2str(scanN); % store scan#
                                bl_MTloc(bl_MTlocRow,5) = {iB}; % store block#

                                bl_MTloc(bl_MTlocRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for moving condition
                                bl_MTloc(bl_MTlocRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for moving condition
                                bl_MTloc(bl_MTlocRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for moving condition
                                bl_MTloc(bl_MTlocRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for moving condition
                                bl_MTloc(bl_MTlocRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''distanceMean''}{1}(iB)}','cell'); % store distance M for moving condition
                                bl_MTloc(bl_MTlocRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for moving condition
                                bl_MTloc(bl_MTlocRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''angleMean''}{1}(iB)}','cell'); % store angle M for moving condition
                                bl_MTloc(bl_MTlocRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''angleSD''}{1}(iB)}','cell'); % store angle SD for moving condition
                                bl_MTloc(bl_MTlocRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for moving condition
                                bl_MTloc(bl_MTlocRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for moving condition
                                bl_MTloc(bl_MTlocRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for moving condition
                                bl_MTloc(bl_MTlocRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for moving condition
                                bl_MTloc(bl_MTlocRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for moving condition
                                bl_MTloc(bl_MTlocRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for moving condition
                                bl_MTloc(bl_MTlocRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for moving condition
                                bl_MTloc(bl_MTlocRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for moving condition
                                bl_MTloc(bl_MTlocRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for moving condition
                                bl_MTloc(bl_MTlocRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for moving condition
                                bl_MTloc(bl_MTlocRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''moving'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for moving condition

                                bl_MTloc(bl_MTlocRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for static condition
                                bl_MTloc(bl_MTlocRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for static condition
                                bl_MTloc(bl_MTlocRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for static condition
                                bl_MTloc(bl_MTlocRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for static condition
                                bl_MTloc(bl_MTlocRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''distanceMean''}{1}(iB)}','cell'); % store distance M for static condition
                                bl_MTloc(bl_MTlocRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for static condition
                                bl_MTloc(bl_MTlocRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''angleMean''}{1}(iB)}','cell'); % store angle M for static condition
                                bl_MTloc(bl_MTlocRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''angleSD''}{1}(iB)}','cell'); % store angle SD for static condition
                                bl_MTloc(bl_MTlocRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for static condition
                                bl_MTloc(bl_MTlocRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for static condition
                                bl_MTloc(bl_MTlocRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for static condition
                                bl_MTloc(bl_MTlocRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for static condition
                                bl_MTloc(bl_MTlocRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for static condition
                                bl_MTloc(bl_MTlocRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for static condition
                                bl_MTloc(bl_MTlocRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for static condition
                                bl_MTloc(bl_MTlocRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for static condition
                                bl_MTloc(bl_MTlocRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for static condition
                                bl_MTloc(bl_MTlocRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for static condition
                                bl_MTloc(bl_MTlocRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''static'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for static condition

                                bl_MTlocRow = bl_MTlocRow+1; % advance counter
                            end

                        elseif AK_whichPattern(data.logFilename{iF,iA},V1_fixlocStr,true)==1 && istable(data.summaryStats{iF,iA}) && data.quality{iF,iA}==1 % is this an V1_fix localizer?; must check before V1 localizer because of overlap in strings
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,V1_fixloc(V1_fixlocRow,5:42))) % check to advance counter
                                V1_fixlocRow = V1_fixlocRow+1; % advance counter
                            elseif any(~cellfun(@isempty,V1_fixloc(V1_fixlocRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(V1_fixlocRow) ' of V1_fixloc']) % message
                                V1_fixlocRow = V1_fixlocRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            V1_fixloc{V1_fixlocRow,1} = subjects{iS}; % store subject name
                            V1_fixloc{V1_fixlocRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            V1_fixloc{V1_fixlocRow,3} = data.setN{iF,iA}; % store set number
                            V1_fixloc{V1_fixlocRow,4} = num2str(scanN); % store scan#

                            V1_fixloc(V1_fixlocRow,5) = table2cell(data.summaryStats{iF,iA}('surr','noTrackTime')); % store no tracking time for fix condition
                            V1_fixloc(V1_fixlocRow,6) = table2cell(data.summaryStats{iF,iA}('surr','fixationTime')); % store fixation time for fix condition
                            V1_fixloc(V1_fixlocRow,7) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedFixationTime')); % store drift-corrected fixation time for fix condition
                            V1_fixloc(V1_fixlocRow,8) = table2cell(data.summaryStats{iF,iA}('surr','saccadeTime')); % store saccade and blink time for fix condition
                            V1_fixloc(V1_fixlocRow,9) = table2cell(data.summaryStats{iF,iA}('surr','distanceMean')); % store distance M for fix condition
                            V1_fixloc(V1_fixlocRow,10) = table2cell(data.summaryStats{iF,iA}('surr','distanceSD')); % store distance SD for fix condition
                            V1_fixloc(V1_fixlocRow,11) = table2cell(data.summaryStats{iF,iA}('surr','angleMean')); % store angle M for fix condition
                            V1_fixloc(V1_fixlocRow,12) = table2cell(data.summaryStats{iF,iA}('surr','angleSD')); % store angle SD for fix condition
                            V1_fixloc(V1_fixlocRow,13) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedDistanceMean')); % store drift-corrected distance M for fix condition
                            V1_fixloc(V1_fixlocRow,14) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedDistanceSD')); % store drift-corrected distance SD for fix condition
                            V1_fixloc(V1_fixlocRow,15) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedAngleMean')); % store drift-corrected angle M for fix condition
                            V1_fixloc(V1_fixlocRow,16) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedAngleSD')); % store drift-corrected angle SD for fix condition
                            V1_fixloc(V1_fixlocRow,17) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeN')); % store number of saccades for fix condition
                            V1_fixloc(V1_fixlocRow,18) = table2cell(data.summaryStats{iF,iA}('surr','bigSaccadeN')); % store number of big saccades for fix condition
                            V1_fixloc(V1_fixlocRow,19) = table2cell(data.summaryStats{iF,iA}('surr','BlinkN')); % store number of saccades for fix condition
                            V1_fixloc(V1_fixlocRow,20) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeAmplitudeMean')); % store saccade amplitude mean for fix condition
                            V1_fixloc(V1_fixlocRow,21) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeAmplitudeSD')); % store saccade amplitude SD for fix condition
                            V1_fixloc(V1_fixlocRow,22) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeVelocityMean')); % store saccade velocity mean for fix condition
                            V1_fixloc(V1_fixlocRow,23) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeVelocitySD')); % store saccade velocity SD for fix condition

                            V1_fixloc(V1_fixlocRow,24) = table2cell(data.summaryStats{iF,iA}('tar','noTrackTime')); % store no tracking time for small condition
                            V1_fixloc(V1_fixlocRow,25) = table2cell(data.summaryStats{iF,iA}('tar','fixationTime')); % store fixation time for small condition
                            V1_fixloc(V1_fixlocRow,26) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedFixationTime')); % store drift-corrected fixation time for small condition
                            V1_fixloc(V1_fixlocRow,27) = table2cell(data.summaryStats{iF,iA}('tar','saccadeTime')); % store saccade and blink time for small condition
                            V1_fixloc(V1_fixlocRow,28) = table2cell(data.summaryStats{iF,iA}('tar','distanceMean')); % store distance M for small condition
                            V1_fixloc(V1_fixlocRow,29) = table2cell(data.summaryStats{iF,iA}('tar','distanceSD')); % store distance SD for small condition
                            V1_fixloc(V1_fixlocRow,30) = table2cell(data.summaryStats{iF,iA}('tar','angleMean')); % store angle M for small condition
                            V1_fixloc(V1_fixlocRow,31) = table2cell(data.summaryStats{iF,iA}('tar','angleSD')); % store angle SD for small condition
                            V1_fixloc(V1_fixlocRow,32) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedDistanceMean')); % store drift-corrected distance M for small condition
                            V1_fixloc(V1_fixlocRow,33) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedDistanceSD')); % store drift-corrected distance SD for small condition
                            V1_fixloc(V1_fixlocRow,34) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedAngleMean')); % store drift-corrected angle M for small condition
                            V1_fixloc(V1_fixlocRow,35) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedAngleSD')); % store drift-corrected angle SD for small condition
                            V1_fixloc(V1_fixlocRow,36) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeN')); % store number of saccades for small condition
                            V1_fixloc(V1_fixlocRow,37) = table2cell(data.summaryStats{iF,iA}('tar','bigSaccadeN')); % store number of big saccades for small condition
                            V1_fixloc(V1_fixlocRow,38) = table2cell(data.summaryStats{iF,iA}('tar','BlinkN')); % store number of saccades for small condition
                            V1_fixloc(V1_fixlocRow,39) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeAmplitudeMean')); % store saccade amplitude mean for small condition
                            V1_fixloc(V1_fixlocRow,40) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeAmplitudeSD')); % store saccade amplitude SD for small condition
                            V1_fixloc(V1_fixlocRow,41) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeVelocityMean')); % store saccade velocity mean for small condition
                            V1_fixloc(V1_fixlocRow,42) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeVelocitySD')); % store saccade velocity SD for small condition

                            V1_fixloc(V1_fixlocRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            V1_fixloc(V1_fixlocRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            V1_fixloc(V1_fixlocRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            V1_fixloc(V1_fixlocRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            V1_fixloc(V1_fixlocRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            V1_fixloc(V1_fixlocRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            V1_fixloc(V1_fixlocRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            V1_fixloc(V1_fixlocRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            V1_fixloc(V1_fixlocRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            V1_fixloc(V1_fixlocRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            V1_fixloc(V1_fixlocRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            V1_fixloc(V1_fixlocRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            V1_fixloc(V1_fixlocRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            V1_fixloc(V1_fixlocRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            V1_fixloc(V1_fixlocRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            V1_fixloc(V1_fixlocRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            V1_fixloc(V1_fixlocRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            V1_fixloc(V1_fixlocRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            V1_fixloc(V1_fixlocRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_V1_fixloc{bl_V1_fixlocRow,1} = subjects{iS}; % store subject name
                                bl_V1_fixloc{bl_V1_fixlocRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_V1_fixloc{bl_V1_fixlocRow,3} = data.setN{iF,iA}; % store set number
                                bl_V1_fixloc{bl_V1_fixlocRow,4} = num2str(scanN); % store scan#
                                bl_V1_fixloc(bl_V1_fixlocRow,5) = {iB}; % store block#

                                bl_V1_fixloc(bl_V1_fixlocRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''distanceMean''}{1}(iB)}','cell'); % store distance M for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''angleMean''}{1}(iB)}','cell'); % store angle M for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''angleSD''}{1}(iB)}','cell'); % store angle SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for fix condition
                                bl_V1_fixloc(bl_V1_fixlocRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for fix condition

                                bl_V1_fixloc(bl_V1_fixlocRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''distanceMean''}{1}(iB)}','cell'); % store distance M for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''angleMean''}{1}(iB)}','cell'); % store angle M for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''angleSD''}{1}(iB)}','cell'); % store angle SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for small condition
                                bl_V1_fixloc(bl_V1_fixlocRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for small condition

                                bl_V1_fixlocRow = bl_V1_fixlocRow+1; % advance counter
                            end

                        elseif AK_whichPattern(data.logFilename{iF,iA},V1locStr,true)==1 && istable(data.summaryStats{iF,iA}) && data.quality{iF,iA}==1 % is this an V1 localizer?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,V1loc(V1locRow,5:42))) % check to advance counter
                                V1locRow = V1locRow+1; % advance counter
                            elseif any(~cellfun(@isempty,V1loc(V1locRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(V1locRow) ' of V1loc']) % message
                                V1locRow = V1locRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            V1loc{V1locRow,1} = subjects{iS}; % store subject name
                            V1loc{V1locRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            V1loc{V1locRow,3} = data.setN{iF,iA}; % store set number
                            V1loc{V1locRow,4} = num2str(scanN); % store scan#

                            V1loc(V1locRow,5) = table2cell(data.summaryStats{iF,iA}('surr','noTrackTime')); % store no tracking time for big condition
                            V1loc(V1locRow,6) = table2cell(data.summaryStats{iF,iA}('surr','fixationTime')); % store fixation time for big condition
                            V1loc(V1locRow,7) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedFixationTime')); % store drift-corrected fixation time for big condition
                            V1loc(V1locRow,8) = table2cell(data.summaryStats{iF,iA}('surr','saccadeTime')); % store saccade and blink time for big condition
                            V1loc(V1locRow,9) = table2cell(data.summaryStats{iF,iA}('surr','distanceMean')); % store distance M for big condition
                            V1loc(V1locRow,10) = table2cell(data.summaryStats{iF,iA}('surr','distanceSD')); % store distance SD for big condition
                            V1loc(V1locRow,11) = table2cell(data.summaryStats{iF,iA}('surr','angleMean')); % store angle M for big condition
                            V1loc(V1locRow,12) = table2cell(data.summaryStats{iF,iA}('surr','angleSD')); % store angle SD for big condition
                            V1loc(V1locRow,13) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedDistanceMean')); % store drift-corrected distance M for big condition
                            V1loc(V1locRow,14) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedDistanceSD')); % store drift-corrected distance SD for big condition
                            V1loc(V1locRow,15) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedAngleMean')); % store drift-corrected angle M for big condition
                            V1loc(V1locRow,16) = table2cell(data.summaryStats{iF,iA}('surr','driftCorrectedAngleSD')); % store drift-corrected angle SD for big condition
                            V1loc(V1locRow,17) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeN')); % store number of saccades for big condition
                            V1loc(V1locRow,18) = table2cell(data.summaryStats{iF,iA}('surr','bigSaccadeN')); % store number of big saccades for big condition
                            V1loc(V1locRow,19) = table2cell(data.summaryStats{iF,iA}('surr','BlinkN')); % store number of saccades for big condition
                            V1loc(V1locRow,20) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeAmplitudeMean')); % store saccade amplitude mean for big condition
                            V1loc(V1locRow,21) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeAmplitudeSD')); % store saccade amplitude SD for big condition
                            V1loc(V1locRow,22) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeVelocityMean')); % store saccade velocity mean for big condition
                            V1loc(V1locRow,23) = table2cell(data.summaryStats{iF,iA}('surr','SaccadeVelocitySD')); % store saccade velocity SD for big condition

                            V1loc(V1locRow,24) = table2cell(data.summaryStats{iF,iA}('tar','noTrackTime')); % store no tracking time for small condition
                            V1loc(V1locRow,25) = table2cell(data.summaryStats{iF,iA}('tar','fixationTime')); % store fixation time for small condition
                            V1loc(V1locRow,26) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedFixationTime')); % store drift-corrected fixation time for small condition
                            V1loc(V1locRow,27) = table2cell(data.summaryStats{iF,iA}('tar','saccadeTime')); % store saccade and blink time for small condition
                            V1loc(V1locRow,28) = table2cell(data.summaryStats{iF,iA}('tar','distanceMean')); % store distance M for small condition
                            V1loc(V1locRow,29) = table2cell(data.summaryStats{iF,iA}('tar','distanceSD')); % store distance SD for small condition
                            V1loc(V1locRow,30) = table2cell(data.summaryStats{iF,iA}('tar','angleMean')); % store angle M for small condition
                            V1loc(V1locRow,31) = table2cell(data.summaryStats{iF,iA}('tar','angleSD')); % store angle SD for small condition
                            V1loc(V1locRow,32) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedDistanceMean')); % store drift-corrected distance M for small condition
                            V1loc(V1locRow,33) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedDistanceSD')); % store drift-corrected distance SD for small condition
                            V1loc(V1locRow,34) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedAngleMean')); % store drift-corrected angle M for small condition
                            V1loc(V1locRow,35) = table2cell(data.summaryStats{iF,iA}('tar','driftCorrectedAngleSD')); % store drift-corrected angle SD for small condition
                            V1loc(V1locRow,36) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeN')); % store number of saccades for small condition
                            V1loc(V1locRow,37) = table2cell(data.summaryStats{iF,iA}('tar','bigSaccadeN')); % store number of big saccades for small condition
                            V1loc(V1locRow,38) = table2cell(data.summaryStats{iF,iA}('tar','BlinkN')); % store number of saccades for small condition
                            V1loc(V1locRow,39) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeAmplitudeMean')); % store saccade amplitude mean for small condition
                            V1loc(V1locRow,40) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeAmplitudeSD')); % store saccade amplitude SD for small condition
                            V1loc(V1locRow,41) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeVelocityMean')); % store saccade velocity mean for small condition
                            V1loc(V1locRow,42) = table2cell(data.summaryStats{iF,iA}('tar','SaccadeVelocitySD')); % store saccade velocity SD for small condition

                            V1loc(V1locRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            V1loc(V1locRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            V1loc(V1locRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            V1loc(V1locRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            V1loc(V1locRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            V1loc(V1locRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            V1loc(V1locRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            V1loc(V1locRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            V1loc(V1locRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            V1loc(V1locRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            V1loc(V1locRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            V1loc(V1locRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            V1loc(V1locRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            V1loc(V1locRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            V1loc(V1locRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            V1loc(V1locRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            V1loc(V1locRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            V1loc(V1locRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            V1loc(V1locRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_V1loc{bl_V1locRow,1} = subjects{iS}; % store subject name
                                bl_V1loc{bl_V1locRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_V1loc{bl_V1locRow,3} = data.setN{iF,iA}; % store set number
                                bl_V1loc{bl_V1locRow,4} = num2str(scanN); % store scan#
                                bl_V1loc(bl_V1locRow,5) = {iB}; % store block#

                                bl_V1loc(bl_V1locRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for big condition
                                bl_V1loc(bl_V1locRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for big condition
                                bl_V1loc(bl_V1locRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for big condition
                                bl_V1loc(bl_V1locRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for big condition
                                bl_V1loc(bl_V1locRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''distanceMean''}{1}(iB)}','cell'); % store distance M for big condition
                                bl_V1loc(bl_V1locRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for big condition
                                bl_V1loc(bl_V1locRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''angleMean''}{1}(iB)}','cell'); % store angle M for big condition
                                bl_V1loc(bl_V1locRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''angleSD''}{1}(iB)}','cell'); % store angle SD for big condition
                                bl_V1loc(bl_V1locRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for big condition
                                bl_V1loc(bl_V1locRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for big condition
                                bl_V1loc(bl_V1locRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for big condition
                                bl_V1loc(bl_V1locRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for big condition
                                bl_V1loc(bl_V1locRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for big condition
                                bl_V1loc(bl_V1locRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for big condition
                                bl_V1loc(bl_V1locRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for big condition
                                bl_V1loc(bl_V1locRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for big condition
                                bl_V1loc(bl_V1locRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for big condition
                                bl_V1loc(bl_V1locRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for big condition
                                bl_V1loc(bl_V1locRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''surr'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for big condition

                                bl_V1loc(bl_V1locRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for small condition
                                bl_V1loc(bl_V1locRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for small condition
                                bl_V1loc(bl_V1locRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for small condition
                                bl_V1loc(bl_V1locRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for small condition
                                bl_V1loc(bl_V1locRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''distanceMean''}{1}(iB)}','cell'); % store distance M for small condition
                                bl_V1loc(bl_V1locRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for small condition
                                bl_V1loc(bl_V1locRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''angleMean''}{1}(iB)}','cell'); % store angle M for small condition
                                bl_V1loc(bl_V1locRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''angleSD''}{1}(iB)}','cell'); % store angle SD for small condition
                                bl_V1loc(bl_V1locRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for small condition
                                bl_V1loc(bl_V1locRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for small condition
                                bl_V1loc(bl_V1locRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for small condition
                                bl_V1loc(bl_V1locRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for small condition
                                bl_V1loc(bl_V1locRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for small condition
                                bl_V1loc(bl_V1locRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for small condition
                                bl_V1loc(bl_V1locRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for small condition
                                bl_V1loc(bl_V1locRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for small condition
                                bl_V1loc(bl_V1locRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for small condition
                                bl_V1loc(bl_V1locRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for small condition
                                bl_V1loc(bl_V1locRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''tar'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for small condition

                                bl_V1locRow = bl_V1locRow+1; % advance counter
                            end

                        elseif AK_whichPattern(data.logFilename{iF,iA},contrastStr,true)==1 && istable(data.summaryStats{iF,iA}) && data.quality{iF,iA}==1 % is this a contrast scan?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,contrast(contrastRow,5:61))) % check to advance counter
                                contrastRow = contrastRow+1; % advance counter
                            elseif any(~cellfun(@isempty,contrast(contrastRow,5:61))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(contrastRow) ' of contrast']) % message
                                contrastRow = contrastRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            contrast{contrastRow,1} = subjects{iS}; % store subject name
                            contrast{contrastRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            contrast{contrastRow,3} = data.setN{iF,iA}; % store set number
                            contrast{contrastRow,4} = num2str(scanN); % store scan#

                            contrast(contrastRow,5) = table2cell(data.summaryStats{iF,iA}('fix','noTrackTime')); % store no tracking time for fix condition
                            contrast(contrastRow,6) = table2cell(data.summaryStats{iF,iA}('fix','fixationTime')); % store fixation time for fix condition
                            contrast(contrastRow,7) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedFixationTime')); % store drift-corrected fixation time for fix condition
                            contrast(contrastRow,8) = table2cell(data.summaryStats{iF,iA}('fix','saccadeTime')); % store saccade and blink time for fix condition
                            contrast(contrastRow,9) = table2cell(data.summaryStats{iF,iA}('fix','distanceMean')); % store distance M for fix condition
                            contrast(contrastRow,10) = table2cell(data.summaryStats{iF,iA}('fix','distanceSD')); % store distance SD for fix condition
                            contrast(contrastRow,11) = table2cell(data.summaryStats{iF,iA}('fix','angleMean')); % store angle M for fix condition
                            contrast(contrastRow,12) = table2cell(data.summaryStats{iF,iA}('fix','angleSD')); % store angle SD for fix condition
                            contrast(contrastRow,13) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedDistanceMean')); % store drift-corrected distance M for fix condition
                            contrast(contrastRow,14) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedDistanceSD')); % store drift-corrected distance SD for fix condition
                            contrast(contrastRow,15) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedAngleMean')); % store drift-corrected angle M for fix condition
                            contrast(contrastRow,16) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedAngleSD')); % store drift-corrected angle SD for fix condition
                            contrast(contrastRow,17) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeN')); % store number of saccades for fix condition
                            contrast(contrastRow,18) = table2cell(data.summaryStats{iF,iA}('fix','bigSaccadeN')); % store number of big saccades for fix condition
                            contrast(contrastRow,19) = table2cell(data.summaryStats{iF,iA}('fix','BlinkN')); % store number of saccades for fix condition
                            contrast(contrastRow,20) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeAmplitudeMean')); % store saccade amplitude mean for fix condition
                            contrast(contrastRow,21) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeAmplitudeSD')); % store saccade amplitude SD for fix condition
                            contrast(contrastRow,22) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeVelocityMean')); % store saccade velocity mean for fix condition
                            contrast(contrastRow,23) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeVelocitySD')); % store saccade velocity SD for fix condition

                            contrast(contrastRow,24) = table2cell(data.summaryStats{iF,iA}('low','noTrackTime')); % store no tracking time for low condition
                            contrast(contrastRow,25) = table2cell(data.summaryStats{iF,iA}('low','fixationTime')); % store fixation time for low condition
                            contrast(contrastRow,26) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedFixationTime')); % store drift-corrected fixation time for low condition
                            contrast(contrastRow,27) = table2cell(data.summaryStats{iF,iA}('low','saccadeTime')); % store saccade and blink time for low condition
                            contrast(contrastRow,28) = table2cell(data.summaryStats{iF,iA}('low','distanceMean')); % store distance M for low condition
                            contrast(contrastRow,29) = table2cell(data.summaryStats{iF,iA}('low','distanceSD')); % store distance SD for low condition
                            contrast(contrastRow,30) = table2cell(data.summaryStats{iF,iA}('low','angleMean')); % store angle M for low condition
                            contrast(contrastRow,31) = table2cell(data.summaryStats{iF,iA}('low','angleSD')); % store angle SD for low condition
                            contrast(contrastRow,32) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedDistanceMean')); % store drift-corrected distance M for low condition
                            contrast(contrastRow,33) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedDistanceSD')); % store drift-corrected distance SD for low condition
                            contrast(contrastRow,34) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedAngleMean')); % store drift-corrected angle M for low condition
                            contrast(contrastRow,35) = table2cell(data.summaryStats{iF,iA}('low','driftCorrectedAngleSD')); % store drift-corrected angle SD for low condition
                            contrast(contrastRow,36) = table2cell(data.summaryStats{iF,iA}('low','SaccadeN')); % store number of saccades for low condition
                            contrast(contrastRow,37) = table2cell(data.summaryStats{iF,iA}('low','bigSaccadeN')); % store number of big saccades for low condition
                            contrast(contrastRow,38) = table2cell(data.summaryStats{iF,iA}('low','BlinkN')); % store number of saccades for low condition
                            contrast(contrastRow,39) = table2cell(data.summaryStats{iF,iA}('low','SaccadeAmplitudeMean')); % store saccade amplitude mean for low condition
                            contrast(contrastRow,40) = table2cell(data.summaryStats{iF,iA}('low','SaccadeAmplitudeSD')); % store saccade amplitude SD for low condition
                            contrast(contrastRow,41) = table2cell(data.summaryStats{iF,iA}('low','SaccadeVelocityMean')); % store saccade velocity mean for low condition
                            contrast(contrastRow,42) = table2cell(data.summaryStats{iF,iA}('low','SaccadeVelocitySD')); % store saccade velocity SD for low condition

                            contrast(contrastRow,43) = table2cell(data.summaryStats{iF,iA}('hi','noTrackTime')); % store no tracking time for hi condition
                            contrast(contrastRow,44) = table2cell(data.summaryStats{iF,iA}('hi','fixationTime')); % store fixation time for hi condition
                            contrast(contrastRow,45) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedFixationTime')); % store drift-corrected fixation time for hi condition
                            contrast(contrastRow,46) = table2cell(data.summaryStats{iF,iA}('hi','saccadeTime')); % store saccade and blink time for hi condition
                            contrast(contrastRow,47) = table2cell(data.summaryStats{iF,iA}('hi','distanceMean')); % store distance M for hi condition
                            contrast(contrastRow,48) = table2cell(data.summaryStats{iF,iA}('hi','distanceSD')); % store distance SD for hi condition
                            contrast(contrastRow,49) = table2cell(data.summaryStats{iF,iA}('hi','angleMean')); % store angle M for hi condition
                            contrast(contrastRow,50) = table2cell(data.summaryStats{iF,iA}('hi','angleSD')); % store angle SD for hi condition
                            contrast(contrastRow,51) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedDistanceMean')); % store drift-corrected distance M for hi condition
                            contrast(contrastRow,52) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedDistanceSD')); % store drift-corrected distance SD for hi condition
                            contrast(contrastRow,53) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedAngleMean')); % store drift-corrected angle M for hi condition
                            contrast(contrastRow,54) = table2cell(data.summaryStats{iF,iA}('hi','driftCorrectedAngleSD')); % store drift-corrected angle SD for hi condition
                            contrast(contrastRow,55) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeN')); % store number of saccades for hi condition
                            contrast(contrastRow,56) = table2cell(data.summaryStats{iF,iA}('hi','bigSaccadeN')); % store number of big saccades for hi condition
                            contrast(contrastRow,57) = table2cell(data.summaryStats{iF,iA}('hi','BlinkN')); % store number of saccades for hi condition
                            contrast(contrastRow,58) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeAmplitudeMean')); % store saccade amplitude mean for hi condition
                            contrast(contrastRow,59) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeAmplitudeSD')); % store saccade amplitude SD for hi condition
                            contrast(contrastRow,60) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeVelocityMean')); % store saccade velocity mean for hi condition
                            contrast(contrastRow,61) = table2cell(data.summaryStats{iF,iA}('hi','SaccadeVelocitySD')); % store saccade velocity SD for hi condition

                            contrast(contrastRow,62) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            contrast(contrastRow,63) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            contrast(contrastRow,64) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            contrast(contrastRow,65) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            contrast(contrastRow,66) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            contrast(contrastRow,67) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            contrast(contrastRow,68) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            contrast(contrastRow,69) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            contrast(contrastRow,70) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            contrast(contrastRow,71) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            contrast(contrastRow,72) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            contrast(contrastRow,73) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            contrast(contrastRow,74) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            contrast(contrastRow,75) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            contrast(contrastRow,76) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            contrast(contrastRow,77) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            contrast(contrastRow,78) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            contrast(contrastRow,79) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            contrast(contrastRow,80) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_contrast{bl_contrastRow,1} = subjects{iS}; % store subject name
                                bl_contrast{bl_contrastRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_contrast{bl_contrastRow,3} = data.setN{iF,iA}; % store set number
                                bl_contrast{bl_contrastRow,4} = num2str(scanN); % store scan#
                                bl_contrast(bl_contrastRow,5) = {iB}; % store block#

                                bl_contrast(bl_contrastRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for fix condition
                                bl_contrast(bl_contrastRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for fix condition
                                bl_contrast(bl_contrastRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for fix condition
                                bl_contrast(bl_contrastRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for fix condition
                                bl_contrast(bl_contrastRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''distanceMean''}{1}(iB)}','cell'); % store distance M for fix condition
                                bl_contrast(bl_contrastRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for fix condition
                                bl_contrast(bl_contrastRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''angleMean''}{1}(iB)}','cell'); % store angle M for fix condition
                                bl_contrast(bl_contrastRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''angleSD''}{1}(iB)}','cell'); % store angle SD for fix condition
                                bl_contrast(bl_contrastRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for fix condition
                                bl_contrast(bl_contrastRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for fix condition
                                bl_contrast(bl_contrastRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for fix condition
                                bl_contrast(bl_contrastRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for fix condition
                                bl_contrast(bl_contrastRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for fix condition
                                bl_contrast(bl_contrastRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for fix condition
                                bl_contrast(bl_contrastRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for fix condition
                                bl_contrast(bl_contrastRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for fix condition
                                bl_contrast(bl_contrastRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for fix condition
                                bl_contrast(bl_contrastRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for fix condition
                                bl_contrast(bl_contrastRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for fix condition

                                bl_contrast(bl_contrastRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for low condition
                                bl_contrast(bl_contrastRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for low condition
                                bl_contrast(bl_contrastRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for low condition
                                bl_contrast(bl_contrastRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for low condition
                                bl_contrast(bl_contrastRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''distanceMean''}{1}(iB)}','cell'); % store distance M for low condition
                                bl_contrast(bl_contrastRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for low condition
                                bl_contrast(bl_contrastRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''angleMean''}{1}(iB)}','cell'); % store angle M for low condition
                                bl_contrast(bl_contrastRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''angleSD''}{1}(iB)}','cell'); % store angle SD for low condition
                                bl_contrast(bl_contrastRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for low condition
                                bl_contrast(bl_contrastRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for low condition
                                bl_contrast(bl_contrastRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for low condition
                                bl_contrast(bl_contrastRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for low condition
                                bl_contrast(bl_contrastRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for low condition
                                bl_contrast(bl_contrastRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for low condition
                                bl_contrast(bl_contrastRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for low condition
                                bl_contrast(bl_contrastRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for low condition
                                bl_contrast(bl_contrastRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for low condition
                                bl_contrast(bl_contrastRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for low condition
                                bl_contrast(bl_contrastRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''low'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for low condition

                                bl_contrast(bl_contrastRow,44) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for hi condition
                                bl_contrast(bl_contrastRow,45) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for hi condition
                                bl_contrast(bl_contrastRow,46) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for hi condition
                                bl_contrast(bl_contrastRow,47) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for hi condition
                                bl_contrast(bl_contrastRow,48) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''distanceMean''}{1}(iB)}','cell'); % store distance M for hi condition
                                bl_contrast(bl_contrastRow,49) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for hi condition
                                bl_contrast(bl_contrastRow,50) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''angleMean''}{1}(iB)}','cell'); % store angle M for hi condition
                                bl_contrast(bl_contrastRow,51) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''angleSD''}{1}(iB)}','cell'); % store angle SD for hi condition
                                bl_contrast(bl_contrastRow,52) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for hi condition
                                bl_contrast(bl_contrastRow,53) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for hi condition
                                bl_contrast(bl_contrastRow,54) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for hi condition
                                bl_contrast(bl_contrastRow,55) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for hi condition
                                bl_contrast(bl_contrastRow,56) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for hi condition
                                bl_contrast(bl_contrastRow,57) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for hi condition
                                bl_contrast(bl_contrastRow,58) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for hi condition
                                bl_contrast(bl_contrastRow,59) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for hi condition
                                bl_contrast(bl_contrastRow,60) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for hi condition
                                bl_contrast(bl_contrastRow,61) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for hi condition
                                bl_contrast(bl_contrastRow,62) = AK_catchEmpty('{data.blockStats{iF,iA}{''hi'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for hi condition

                                bl_contrastRow = bl_contrastRow+1; % advance counter
                            end
                        elseif AK_whichPattern(data.logFilename{iF,iA},supStr,true)==1 && istable(data.summaryStats{iF,iA}) && data.quality{iF,iA}==1 % is this a suppression scan?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,sup(supRow,5:42))) % check to advance counter
                                supRow = supRow+1; % advance counter
                            elseif any(~cellfun(@isempty,sup(supRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(supRow) ' of sup']) % message
                                supRow = supRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            sup{supRow,1} = subjects{iS}; % store subject name
                            sup{supRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            sup{supRow,3} = data.setN{iF,iA}; % store set number
                            sup{supRow,4} = num2str(scanN); % store scan#

                            sup(supRow,5) = table2cell(data.summaryStats{iF,iA}('small','noTrackTime')); % store no tracking time for small condition
                            sup(supRow,6) = table2cell(data.summaryStats{iF,iA}('small','fixationTime')); % store fixation time for small condition
                            sup(supRow,7) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedFixationTime')); % store drift-corrected fixation time for small condition
                            sup(supRow,8) = table2cell(data.summaryStats{iF,iA}('small','saccadeTime')); % store saccade and blink time for small condition
                            sup(supRow,9) = table2cell(data.summaryStats{iF,iA}('small','distanceMean')); % store distance M for small condition
                            sup(supRow,10) = table2cell(data.summaryStats{iF,iA}('small','distanceSD')); % store distance SD for small condition
                            sup(supRow,11) = table2cell(data.summaryStats{iF,iA}('small','angleMean')); % store angle M for small condition
                            sup(supRow,12) = table2cell(data.summaryStats{iF,iA}('small','angleSD')); % store angle SD for small condition
                            sup(supRow,13) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedDistanceMean')); % store drift-corrected distance M for small condition
                            sup(supRow,14) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedDistanceSD')); % store drift-corrected distance SD for small condition
                            sup(supRow,15) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedAngleMean')); % store drift-corrected angle M for small condition
                            sup(supRow,16) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedAngleSD')); % store drift-corrected angle SD for small condition
                            sup(supRow,17) = table2cell(data.summaryStats{iF,iA}('small','SaccadeN')); % store number of saccades for small condition
                            sup(supRow,18) = table2cell(data.summaryStats{iF,iA}('small','bigSaccadeN')); % store number of big saccades for small condition
                            sup(supRow,19) = table2cell(data.summaryStats{iF,iA}('small','BlinkN')); % store number of saccades for small condition
                            sup(supRow,20) = table2cell(data.summaryStats{iF,iA}('small','SaccadeAmplitudeMean')); % store saccade amplitude mean for small condition
                            sup(supRow,21) = table2cell(data.summaryStats{iF,iA}('small','SaccadeAmplitudeSD')); % store saccade amplitude SD for small condition
                            sup(supRow,22) = table2cell(data.summaryStats{iF,iA}('small','SaccadeVelocityMean')); % store saccade velocity mean for small condition
                            sup(supRow,23) = table2cell(data.summaryStats{iF,iA}('small','SaccadeVelocitySD')); % store saccade velocity SD for small condition

                            sup(supRow,24) = table2cell(data.summaryStats{iF,iA}('big','noTrackTime')); % store no tracking time for big condition
                            sup(supRow,25) = table2cell(data.summaryStats{iF,iA}('big','fixationTime')); % store fixation time for big condition
                            sup(supRow,26) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedFixationTime')); % store drift-corrected fixation time for big condition
                            sup(supRow,27) = table2cell(data.summaryStats{iF,iA}('big','saccadeTime')); % store saccade and blink time for big condition
                            sup(supRow,28) = table2cell(data.summaryStats{iF,iA}('big','distanceMean')); % store distance M for big condition
                            sup(supRow,29) = table2cell(data.summaryStats{iF,iA}('big','distanceSD')); % store distance SD for big condition
                            sup(supRow,30) = table2cell(data.summaryStats{iF,iA}('big','angleMean')); % store angle M for big condition
                            sup(supRow,31) = table2cell(data.summaryStats{iF,iA}('big','angleSD')); % store angle SD for big condition
                            sup(supRow,32) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedDistanceMean')); % store drift-corrected distance M for big condition
                            sup(supRow,33) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedDistanceSD')); % store drift-corrected distance SD for big condition
                            sup(supRow,34) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedAngleMean')); % store drift-corrected angle M for big condition
                            sup(supRow,35) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedAngleSD')); % store drift-corrected angle SD for big condition
                            sup(supRow,36) = table2cell(data.summaryStats{iF,iA}('big','SaccadeN')); % store number of saccades for big condition
                            sup(supRow,37) = table2cell(data.summaryStats{iF,iA}('big','bigSaccadeN')); % store number of big saccades for big condition
                            sup(supRow,38) = table2cell(data.summaryStats{iF,iA}('big','BlinkN')); % store number of saccades for big condition
                            sup(supRow,39) = table2cell(data.summaryStats{iF,iA}('big','SaccadeAmplitudeMean')); % store saccade amplitude mean for big condition
                            sup(supRow,40) = table2cell(data.summaryStats{iF,iA}('big','SaccadeAmplitudeSD')); % store saccade amplitude SD for big condition
                            sup(supRow,41) = table2cell(data.summaryStats{iF,iA}('big','SaccadeVelocityMean')); % store saccade velocity mean for big condition
                            sup(supRow,42) = table2cell(data.summaryStats{iF,iA}('big','SaccadeVelocitySD')); % store saccade velocity SD for big condition

                            sup(supRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            sup(supRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            sup(supRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            sup(supRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            sup(supRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            sup(supRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            sup(supRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            sup(supRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            sup(supRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            sup(supRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            sup(supRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            sup(supRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            sup(supRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            sup(supRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            sup(supRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            sup(supRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            sup(supRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            sup(supRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            sup(supRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_sup{bl_supRow,1} = subjects{iS}; % store subject name
                                bl_sup{bl_supRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_sup{bl_supRow,3} = data.setN{iF,iA}; % store set number
                                bl_sup{bl_supRow,4} = num2str(scanN); % store scan#
                                bl_sup(bl_supRow,5) = {iB}; % store block#

                                bl_sup(bl_supRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for small condition
                                bl_sup(bl_supRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for small condition
                                bl_sup(bl_supRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for small condition
                                bl_sup(bl_supRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for small condition
                                bl_sup(bl_supRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''distanceMean''}{1}(iB)}','cell'); % store distance M for small condition
                                bl_sup(bl_supRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for small condition
                                bl_sup(bl_supRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''angleMean''}{1}(iB)}','cell'); % store angle M for small condition
                                bl_sup(bl_supRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''angleSD''}{1}(iB)}','cell'); % store angle SD for small condition
                                bl_sup(bl_supRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for small condition
                                bl_sup(bl_supRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for small condition
                                bl_sup(bl_supRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for small condition
                                bl_sup(bl_supRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for small condition
                                bl_sup(bl_supRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for small condition
                                bl_sup(bl_supRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for small condition
                                bl_sup(bl_supRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for small condition
                                bl_sup(bl_supRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for small condition
                                bl_sup(bl_supRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for small condition
                                bl_sup(bl_supRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for small condition
                                bl_sup(bl_supRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for small condition

                                bl_sup(bl_supRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for big condition
                                bl_sup(bl_supRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for big condition
                                bl_sup(bl_supRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for big condition
                                bl_sup(bl_supRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for big condition
                                bl_sup(bl_supRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''distanceMean''}{1}(iB)}','cell'); % store distance M for big condition
                                bl_sup(bl_supRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for big condition
                                bl_sup(bl_supRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''angleMean''}{1}(iB)}','cell'); % store angle M for big condition
                                bl_sup(bl_supRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''angleSD''}{1}(iB)}','cell'); % store angle SD for big condition
                                bl_sup(bl_supRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for big condition
                                bl_sup(bl_supRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for big condition
                                bl_sup(bl_supRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for big condition
                                bl_sup(bl_supRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for big condition
                                bl_sup(bl_supRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for big condition
                                bl_sup(bl_supRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for big condition
                                bl_sup(bl_supRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for big condition
                                bl_sup(bl_supRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for big condition
                                bl_sup(bl_supRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for big condition
                                bl_sup(bl_supRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for big condition
                                bl_sup(bl_supRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for big condition

                                bl_supRow = bl_supRow+1; % advance counter
                            end

                        elseif AK_whichPattern(data.logFilename{iF,iA},sumStr,true)==1 && istable(data.summaryStats{iF,iA}) && data.quality{iF,iA}==1 % is this a summation scan?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,sum(sumRow,5:42))) % check to advance counter
                                sumRow = sumRow+1; % advance counter
                            elseif any(~cellfun(@isempty,sum(sumRow,5:42))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(sumRow) ' of sum']) % message
                                sumRow = sumRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            sum{sumRow,1} = subjects{iS}; % store subject name
                            sum{sumRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            sum{sumRow,3} = data.setN{iF,iA}; % store set number
                            sum{sumRow,4} = num2str(scanN); % store scan#

                            sum(sumRow,5) = table2cell(data.summaryStats{iF,iA}('small','noTrackTime')); % store no tracking time for small condition
                            sum(sumRow,6) = table2cell(data.summaryStats{iF,iA}('small','fixationTime')); % store fixation time for small condition
                            sum(sumRow,7) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedFixationTime')); % store drift-corrected fixation time for small condition
                            sum(sumRow,8) = table2cell(data.summaryStats{iF,iA}('small','saccadeTime')); % store saccade and blink time for small condition
                            sum(sumRow,9) = table2cell(data.summaryStats{iF,iA}('small','distanceMean')); % store distance M for small condition
                            sum(sumRow,10) = table2cell(data.summaryStats{iF,iA}('small','distanceSD')); % store distance SD for small condition
                            sum(sumRow,11) = table2cell(data.summaryStats{iF,iA}('small','angleMean')); % store angle M for small condition
                            sum(sumRow,12) = table2cell(data.summaryStats{iF,iA}('small','angleSD')); % store angle SD for small condition
                            sum(sumRow,13) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedDistanceMean')); % store drift-corrected distance M for small condition
                            sum(sumRow,14) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedDistanceSD')); % store drift-corrected distance SD for small condition
                            sum(sumRow,15) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedAngleMean')); % store drift-corrected angle M for small condition
                            sum(sumRow,16) = table2cell(data.summaryStats{iF,iA}('small','driftCorrectedAngleSD')); % store drift-corrected angle SD for small condition
                            sum(sumRow,17) = table2cell(data.summaryStats{iF,iA}('small','SaccadeN')); % store number of saccades for small condition
                            sum(sumRow,18) = table2cell(data.summaryStats{iF,iA}('small','bigSaccadeN')); % store number of big saccades for small condition
                            sum(sumRow,19) = table2cell(data.summaryStats{iF,iA}('small','BlinkN')); % store number of saccades for small condition
                            sum(sumRow,20) = table2cell(data.summaryStats{iF,iA}('small','SaccadeAmplitudeMean')); % store saccade amplitude mean for small condition
                            sum(sumRow,21) = table2cell(data.summaryStats{iF,iA}('small','SaccadeAmplitudeSD')); % store saccade amplitude SD for small condition
                            sum(sumRow,22) = table2cell(data.summaryStats{iF,iA}('small','SaccadeVelocityMean')); % store saccade velocity mean for small condition
                            sum(sumRow,23) = table2cell(data.summaryStats{iF,iA}('small','SaccadeVelocitySD')); % store saccade velocity SD for small condition

                            sum(sumRow,24) = table2cell(data.summaryStats{iF,iA}('big','noTrackTime')); % store no tracking time for big condition
                            sum(sumRow,25) = table2cell(data.summaryStats{iF,iA}('big','fixationTime')); % store fixation time for big condition
                            sum(sumRow,26) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedFixationTime')); % store drift-corrected fixation time for big condition
                            sum(sumRow,27) = table2cell(data.summaryStats{iF,iA}('big','saccadeTime')); % store saccade and blink time for big condition
                            sum(sumRow,28) = table2cell(data.summaryStats{iF,iA}('big','distanceMean')); % store distance M for big condition
                            sum(sumRow,29) = table2cell(data.summaryStats{iF,iA}('big','distanceSD')); % store distance SD for big condition
                            sum(sumRow,30) = table2cell(data.summaryStats{iF,iA}('big','angleMean')); % store angle M for big condition
                            sum(sumRow,31) = table2cell(data.summaryStats{iF,iA}('big','angleSD')); % store angle SD for big condition
                            sum(sumRow,32) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedDistanceMean')); % store drift-corrected distance M for big condition
                            sum(sumRow,33) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedDistanceSD')); % store drift-corrected distance SD for big condition
                            sum(sumRow,34) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedAngleMean')); % store drift-corrected angle M for big condition
                            sum(sumRow,35) = table2cell(data.summaryStats{iF,iA}('big','driftCorrectedAngleSD')); % store drift-corrected angle SD for big condition
                            sum(sumRow,36) = table2cell(data.summaryStats{iF,iA}('big','SaccadeN')); % store number of saccades for big condition
                            sum(sumRow,37) = table2cell(data.summaryStats{iF,iA}('big','bigSaccadeN')); % store number of big saccades for big condition
                            sum(sumRow,38) = table2cell(data.summaryStats{iF,iA}('big','BlinkN')); % store number of saccades for big condition
                            sum(sumRow,39) = table2cell(data.summaryStats{iF,iA}('big','SaccadeAmplitudeMean')); % store saccade amplitude mean for big condition
                            sum(sumRow,40) = table2cell(data.summaryStats{iF,iA}('big','SaccadeAmplitudeSD')); % store saccade amplitude SD for big condition
                            sum(sumRow,41) = table2cell(data.summaryStats{iF,iA}('big','SaccadeVelocityMean')); % store saccade velocity mean for big condition
                            sum(sumRow,42) = table2cell(data.summaryStats{iF,iA}('big','SaccadeVelocitySD')); % store saccade velocity SD for big condition

                            sum(sumRow,43) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            sum(sumRow,44) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            sum(sumRow,45) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            sum(sumRow,46) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            sum(sumRow,47) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            sum(sumRow,48) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            sum(sumRow,49) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            sum(sumRow,50) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            sum(sumRow,51) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            sum(sumRow,52) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            sum(sumRow,53) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            sum(sumRow,54) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            sum(sumRow,55) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            sum(sumRow,56) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            sum(sumRow,57) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            sum(sumRow,58) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            sum(sumRow,59) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            sum(sumRow,60) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            sum(sumRow,61) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_sum{bl_sumRow,1} = subjects{iS}; % store subject name
                                bl_sum{bl_sumRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_sum{bl_sumRow,3} = data.setN{iF,iA}; % store set number
                                bl_sum{bl_sumRow,4} = num2str(scanN); % store scan#
                                bl_sum(bl_sumRow,5) = {iB}; % store block#

                                bl_sum(bl_sumRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for small condition
                                bl_sum(bl_sumRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for small condition
                                bl_sum(bl_sumRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for small condition
                                bl_sum(bl_sumRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for small condition
                                bl_sum(bl_sumRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''distanceMean''}{1}(iB)}','cell'); % store distance M for small condition
                                bl_sum(bl_sumRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for small condition
                                bl_sum(bl_sumRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''angleMean''}{1}(iB)}','cell'); % store angle M for small condition
                                bl_sum(bl_sumRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''angleSD''}{1}(iB)}','cell'); % store angle SD for small condition
                                bl_sum(bl_sumRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for small condition
                                bl_sum(bl_sumRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for small condition
                                bl_sum(bl_sumRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for small condition
                                bl_sum(bl_sumRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for small condition
                                bl_sum(bl_sumRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for small condition
                                bl_sum(bl_sumRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for small condition
                                bl_sum(bl_sumRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for small condition
                                bl_sum(bl_sumRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for small condition
                                bl_sum(bl_sumRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for small condition
                                bl_sum(bl_sumRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for small condition
                                bl_sum(bl_sumRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''small'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for small condition

                                bl_sum(bl_sumRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for big condition
                                bl_sum(bl_sumRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for big condition
                                bl_sum(bl_sumRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for big condition
                                bl_sum(bl_sumRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for big condition
                                bl_sum(bl_sumRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''distanceMean''}{1}(iB)}','cell'); % store distance M for big condition
                                bl_sum(bl_sumRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for big condition
                                bl_sum(bl_sumRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''angleMean''}{1}(iB)}','cell'); % store angle M for big condition
                                bl_sum(bl_sumRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''angleSD''}{1}(iB)}','cell'); % store angle SD for big condition
                                bl_sum(bl_sumRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for big condition
                                bl_sum(bl_sumRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for big condition
                                bl_sum(bl_sumRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for big condition
                                bl_sum(bl_sumRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for big condition
                                bl_sum(bl_sumRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for big condition
                                bl_sum(bl_sumRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for big condition
                                bl_sum(bl_sumRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for big condition
                                bl_sum(bl_sumRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for big condition
                                bl_sum(bl_sumRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for big condition
                                bl_sum(bl_sumRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for big condition
                                bl_sum(bl_sumRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''big'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for big condition

                                bl_sumRow = bl_sumRow+1; % advance counter
                            end
                            
                         elseif AK_whichPattern(data.logFilename{iF,iA},ftapStr,true)==1 && istable(data.summaryStats{iF,iA}) && data.quality{iF,iA}==1 % is this a ftap scan?
                            % summary stats:
                            % checks to advance row counter and avoid overwriting
                            if all(~cellfun(@isempty,ftap(ftapRow,5:61))) % check to advance counter
                                ftapRow = ftapRow+1; % advance counter
                            elseif any(~cellfun(@isempty,ftap(ftapRow,5:61))) % check that data is not being overwritten
                                disp(['Avoided overwriting for asc file: ' data.ascFilename{iF,iA} ' at row ' num2str(ftapRow) ' of ftap']) % message
                                ftapRow = ftapRow+1; % advance counter
                            end

                            % store stats about eye tracking data
                            ftap{ftapRow,1} = subjects{iS}; % store subject name
                            ftap{ftapRow,2} = data.fMRIsession{iF}; % store fMRI session name
                            ftap{ftapRow,3} = data.setN{iF,iA}; % store set number
                            ftap{ftapRow,4} = num2str(scanN); % store scan#

                            ftap(ftapRow,5) = table2cell(data.summaryStats{iF,iA}('fix','noTrackTime')); % store no tracking time for fix condition
                            ftap(ftapRow,6) = table2cell(data.summaryStats{iF,iA}('fix','fixationTime')); % store fixation time for fix condition
                            ftap(ftapRow,7) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedFixationTime')); % store drift-corrected fixation time for fix condition
                            ftap(ftapRow,8) = table2cell(data.summaryStats{iF,iA}('fix','saccadeTime')); % store saccade and blink time for fix condition
                            ftap(ftapRow,9) = table2cell(data.summaryStats{iF,iA}('fix','distanceMean')); % store distance M for fix condition
                            ftap(ftapRow,10) = table2cell(data.summaryStats{iF,iA}('fix','distanceSD')); % store distance SD for fix condition
                            ftap(ftapRow,11) = table2cell(data.summaryStats{iF,iA}('fix','angleMean')); % store angle M for fix condition
                            ftap(ftapRow,12) = table2cell(data.summaryStats{iF,iA}('fix','angleSD')); % store angle SD for fix condition
                            ftap(ftapRow,13) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedDistanceMean')); % store drift-corrected distance M for fix condition
                            ftap(ftapRow,14) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedDistanceSD')); % store drift-corrected distance SD for fix condition
                            ftap(ftapRow,15) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedAngleMean')); % store drift-corrected angle M for fix condition
                            ftap(ftapRow,16) = table2cell(data.summaryStats{iF,iA}('fix','driftCorrectedAngleSD')); % store drift-corrected angle SD for fix condition
                            ftap(ftapRow,17) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeN')); % store number of saccades for fix condition
                            ftap(ftapRow,18) = table2cell(data.summaryStats{iF,iA}('fix','bigSaccadeN')); % store number of big saccades for fix condition
                            ftap(ftapRow,19) = table2cell(data.summaryStats{iF,iA}('fix','BlinkN')); % store number of saccades for fix condition
                            ftap(ftapRow,20) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeAmplitudeMean')); % store saccade amplitude mean for fix condition
                            ftap(ftapRow,21) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeAmplitudeSD')); % store saccade amplitude SD for fix condition
                            ftap(ftapRow,22) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeVelocityMean')); % store saccade velocity mean for fix condition
                            ftap(ftapRow,23) = table2cell(data.summaryStats{iF,iA}('fix','SaccadeVelocitySD')); % store saccade velocity SD for fix condition

                            ftap(ftapRow,24) = table2cell(data.summaryStats{iF,iA}('standard','noTrackTime')); % store no tracking time for predictable condition
                            ftap(ftapRow,25) = table2cell(data.summaryStats{iF,iA}('standard','fixationTime')); % store fixation time for predictable condition
                            ftap(ftapRow,26) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedFixationTime')); % store drift-corrected fixation time for predictable condition
                            ftap(ftapRow,27) = table2cell(data.summaryStats{iF,iA}('standard','saccadeTime')); % store saccade and blink time for predictable condition
                            ftap(ftapRow,28) = table2cell(data.summaryStats{iF,iA}('standard','distanceMean')); % store distance M for predictable condition
                            ftap(ftapRow,29) = table2cell(data.summaryStats{iF,iA}('standard','distanceSD')); % store distance SD for predictable condition
                            ftap(ftapRow,30) = table2cell(data.summaryStats{iF,iA}('standard','angleMean')); % store angle M for predictable condition
                            ftap(ftapRow,31) = table2cell(data.summaryStats{iF,iA}('standard','angleSD')); % store angle SD for predictable condition
                            ftap(ftapRow,32) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedDistanceMean')); % store drift-corrected distance M for predictable condition
                            ftap(ftapRow,33) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedDistanceSD')); % store drift-corrected distance SD for predictable condition
                            ftap(ftapRow,34) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedAngleMean')); % store drift-corrected angle M for predictable condition
                            ftap(ftapRow,35) = table2cell(data.summaryStats{iF,iA}('standard','driftCorrectedAngleSD')); % store drift-corrected angle SD for predictable condition
                            ftap(ftapRow,36) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeN')); % store number of saccades for predictable condition
                            ftap(ftapRow,37) = table2cell(data.summaryStats{iF,iA}('standard','bigSaccadeN')); % store number of big saccades for predictable condition
                            ftap(ftapRow,38) = table2cell(data.summaryStats{iF,iA}('standard','BlinkN')); % store number of saccades for predictable condition
                            ftap(ftapRow,39) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeAmplitudeMean')); % store saccade amplitude mean for predictable condition
                            ftap(ftapRow,40) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeAmplitudeSD')); % store saccade amplitude SD for predictable condition
                            ftap(ftapRow,41) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeVelocityMean')); % store saccade velocity mean for predictable condition
                            ftap(ftapRow,42) = table2cell(data.summaryStats{iF,iA}('standard','SaccadeVelocitySD')); % store saccade velocity SD for predictable condition

                            ftap(ftapRow,43) = table2cell(data.summaryStats{iF,iA}('variable','noTrackTime')); % store no tracking time for variable condition
                            ftap(ftapRow,44) = table2cell(data.summaryStats{iF,iA}('variable','fixationTime')); % store fixation time for variable condition
                            ftap(ftapRow,45) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedFixationTime')); % store drift-corrected fixation time for variable condition
                            ftap(ftapRow,46) = table2cell(data.summaryStats{iF,iA}('variable','saccadeTime')); % store saccade and blink time for variable condition
                            ftap(ftapRow,47) = table2cell(data.summaryStats{iF,iA}('variable','distanceMean')); % store distance M for variable condition
                            ftap(ftapRow,48) = table2cell(data.summaryStats{iF,iA}('variable','distanceSD')); % store distance SD for variable condition
                            ftap(ftapRow,49) = table2cell(data.summaryStats{iF,iA}('variable','angleMean')); % store angle M for variable condition
                            ftap(ftapRow,50) = table2cell(data.summaryStats{iF,iA}('variable','angleSD')); % store angle SD for variable condition
                            ftap(ftapRow,51) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedDistanceMean')); % store drift-corrected distance M for variable condition
                            ftap(ftapRow,52) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedDistanceSD')); % store drift-corrected distance SD for variable condition
                            ftap(ftapRow,53) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedAngleMean')); % store drift-corrected angle M for variable condition
                            ftap(ftapRow,54) = table2cell(data.summaryStats{iF,iA}('variable','driftCorrectedAngleSD')); % store drift-corrected angle SD for variable condition
                            ftap(ftapRow,55) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeN')); % store number of saccades for variable condition
                            ftap(ftapRow,56) = table2cell(data.summaryStats{iF,iA}('variable','bigSaccadeN')); % store number of big saccades for variable condition
                            ftap(ftapRow,57) = table2cell(data.summaryStats{iF,iA}('variable','BlinkN')); % store number of saccades for variable condition
                            ftap(ftapRow,58) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeAmplitudeMean')); % store saccade amplitude mean for variable condition
                            ftap(ftapRow,59) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeAmplitudeSD')); % store saccade amplitude SD for variable condition
                            ftap(ftapRow,60) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeVelocityMean')); % store saccade velocity mean for variable condition
                            ftap(ftapRow,61) = table2cell(data.summaryStats{iF,iA}('variable','SaccadeVelocitySD')); % store saccade velocity SD for variable condition

                            ftap(ftapRow,62) = table2cell(data.summaryStats{iF,iA}('all','noTrackTime')); % store no tracking time for all conditions
                            ftap(ftapRow,63) = table2cell(data.summaryStats{iF,iA}('all','fixationTime')); % store fixation time for all conditions
                            ftap(ftapRow,64) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedFixationTime')); % store drift-corrected fixation time for all conditions
                            ftap(ftapRow,65) = table2cell(data.summaryStats{iF,iA}('all','saccadeTime')); % store saccade and blink time for all conditions
                            ftap(ftapRow,66) = table2cell(data.summaryStats{iF,iA}('all','distanceMean')); % store distance M for all conditions
                            ftap(ftapRow,67) = table2cell(data.summaryStats{iF,iA}('all','distanceSD')); % store distance SD for all conditions
                            ftap(ftapRow,68) = table2cell(data.summaryStats{iF,iA}('all','angleMean')); % store angle M for all conditions
                            ftap(ftapRow,69) = table2cell(data.summaryStats{iF,iA}('all','angleSD')); % store angle SD for all conditions
                            ftap(ftapRow,70) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceMean')); % store drift-corrected distance M for all conditions
                            ftap(ftapRow,71) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedDistanceSD')); % store drift-corrected distance SD for all conditions
                            ftap(ftapRow,72) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleMean')); % store drift-corrected angle M for all conditions
                            ftap(ftapRow,73) = table2cell(data.summaryStats{iF,iA}('all','driftCorrectedAngleSD')); % store drift-corrected angle SD for all conditions
                            ftap(ftapRow,74) = table2cell(data.summaryStats{iF,iA}('all','SaccadeN')); % store number of saccades for all conditions
                            ftap(ftapRow,75) = table2cell(data.summaryStats{iF,iA}('all','bigSaccadeN')); % store number of big saccades for all conditions
                            ftap(ftapRow,76) = table2cell(data.summaryStats{iF,iA}('all','BlinkN')); % store number of saccades for all conditions
                            ftap(ftapRow,77) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeMean')); % store saccade amplitude mean for all conditions
                            ftap(ftapRow,78) = table2cell(data.summaryStats{iF,iA}('all','SaccadeAmplitudeSD')); % store saccade amplitude SD for all conditions
                            ftap(ftapRow,79) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocityMean')); % store saccade velocity mean for all conditions
                            ftap(ftapRow,80) = table2cell(data.summaryStats{iF,iA}('all','SaccadeVelocitySD')); % store saccade velocity SD for all conditions
                            
                            % block stats:                        
                            % cycle through block number filling in cell array 
                            for iB = 1:max(data.conditionBlockCount{iF,iA})
                                bl_ftap{bl_ftapRow,1} = subjects{iS}; % store subject name
                                bl_ftap{bl_ftapRow,2} = data.fMRIsession{iF}; % store fMRI session name
                                bl_ftap{bl_ftapRow,3} = data.setN{iF,iA}; % store set number
                                bl_ftap{bl_ftapRow,4} = num2str(scanN); % store scan#
                                bl_ftap(bl_ftapRow,5) = {iB}; % store block#

                                bl_ftap(bl_ftapRow,6) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for fix condition
                                bl_ftap(bl_ftapRow,7) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for fix condition
                                bl_ftap(bl_ftapRow,8) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for fix condition
                                bl_ftap(bl_ftapRow,9) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for fix condition
                                bl_ftap(bl_ftapRow,10) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''distanceMean''}{1}(iB)}','cell'); % store distance M for fix condition
                                bl_ftap(bl_ftapRow,11) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for fix condition
                                bl_ftap(bl_ftapRow,12) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''angleMean''}{1}(iB)}','cell'); % store angle M for fix condition
                                bl_ftap(bl_ftapRow,13) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''angleSD''}{1}(iB)}','cell'); % store angle SD for fix condition
                                bl_ftap(bl_ftapRow,14) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for fix condition
                                bl_ftap(bl_ftapRow,15) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for fix condition
                                bl_ftap(bl_ftapRow,16) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for fix condition
                                bl_ftap(bl_ftapRow,17) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for fix condition
                                bl_ftap(bl_ftapRow,18) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for fix condition
                                bl_ftap(bl_ftapRow,19) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for fix condition
                                bl_ftap(bl_ftapRow,20) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for fix condition
                                bl_ftap(bl_ftapRow,21) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for fix condition
                                bl_ftap(bl_ftapRow,22) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for fix condition
                                bl_ftap(bl_ftapRow,23) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for fix condition
                                bl_ftap(bl_ftapRow,24) = AK_catchEmpty('{data.blockStats{iF,iA}{''fix'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for fix condition

                                bl_ftap(bl_ftapRow,25) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for predictable condition
                                bl_ftap(bl_ftapRow,26) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for predictable condition
                                bl_ftap(bl_ftapRow,27) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for predictable condition
                                bl_ftap(bl_ftapRow,28) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for predictable condition
                                bl_ftap(bl_ftapRow,29) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''distanceMean''}{1}(iB)}','cell'); % store distance M for predictable condition
                                bl_ftap(bl_ftapRow,30) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for predictable condition
                                bl_ftap(bl_ftapRow,31) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''angleMean''}{1}(iB)}','cell'); % store angle M for predictable condition
                                bl_ftap(bl_ftapRow,32) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''angleSD''}{1}(iB)}','cell'); % store angle SD for predictable condition
                                bl_ftap(bl_ftapRow,33) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for predictable condition
                                bl_ftap(bl_ftapRow,34) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for predictable condition
                                bl_ftap(bl_ftapRow,35) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for predictable condition
                                bl_ftap(bl_ftapRow,36) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for predictable condition
                                bl_ftap(bl_ftapRow,37) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for predictable condition
                                bl_ftap(bl_ftapRow,38) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for predictable condition
                                bl_ftap(bl_ftapRow,39) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for predictable condition
                                bl_ftap(bl_ftapRow,40) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for predictable condition
                                bl_ftap(bl_ftapRow,41) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for predictable condition
                                bl_ftap(bl_ftapRow,42) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for predictable condition
                                bl_ftap(bl_ftapRow,43) = AK_catchEmpty('{data.blockStats{iF,iA}{''standard'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for predictable condition

                                bl_ftap(bl_ftapRow,44) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''noTrackTime''}{1}(iB)}','cell'); % store no tracking time for variable condition
                                bl_ftap(bl_ftapRow,45) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''fixationTime''}{1}(iB)}','cell'); % store fixation time for variable condition
                                bl_ftap(bl_ftapRow,46) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedFixationTime''}{1}(iB)}','cell'); % store drift-corrected fixation time for variable condition
                                bl_ftap(bl_ftapRow,47) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''saccadeTime''}{1}(iB)}','cell'); % store saccade and blink time for variable condition
                                bl_ftap(bl_ftapRow,48) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''distanceMean''}{1}(iB)}','cell'); % store distance M for variable condition
                                bl_ftap(bl_ftapRow,49) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''distanceSD''}{1}(iB)}','cell'); % store distance SD for variable condition
                                bl_ftap(bl_ftapRow,50) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''angleMean''}{1}(iB)}','cell'); % store angle M for variable condition
                                bl_ftap(bl_ftapRow,51) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''angleSD''}{1}(iB)}','cell'); % store angle SD for variable condition
                                bl_ftap(bl_ftapRow,52) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedDistanceMean''}{1}(iB)}','cell'); % store drift-corrected distance M for variable condition
                                bl_ftap(bl_ftapRow,53) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedDistanceSD''}{1}(iB)}','cell'); % store drift-corrected distance SD for variable condition
                                bl_ftap(bl_ftapRow,54) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedAngleMean''}{1}(iB)}','cell'); % store drift-corrected angle M for variable condition
                                bl_ftap(bl_ftapRow,55) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''driftCorrectedAngleSD''}{1}(iB)}','cell'); % store drift-corrected angle SD for variable condition
                                bl_ftap(bl_ftapRow,56) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeN''}{1}(iB)}','cell'); % store number of saccades for variable condition
                                bl_ftap(bl_ftapRow,57) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''bigSaccadeN''}{1}(iB)}','cell'); % store number of big saccades for variable condition
                                bl_ftap(bl_ftapRow,58) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''BlinkN''}{1}(iB)}','cell'); % store number of blinks for variable condition
                                bl_ftap(bl_ftapRow,59) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeAmplitudeMean''}{1}(iB)}','cell'); % store saccade amplitude mean for variable condition
                                bl_ftap(bl_ftapRow,60) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeAmplitudeSD''}{1}(iB)}','cell'); % store saccade amplitude SD for variable condition
                                bl_ftap(bl_ftapRow,61) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeVelocityMean''}{1}(iB)}','cell'); % store saccade velocity mean for variable condition
                                bl_ftap(bl_ftapRow,62) = AK_catchEmpty('{data.blockStats{iF,iA}{''variable'',''SaccadeVelocitySD''}{1}(iB)}','cell'); % store saccade velocity SD for variable condition

                                bl_ftapRow = bl_ftapRow+1; % advance counter
                            end
                        end
                end
            end

            iF = iF+1; % increase counter
        end
    end
end

% trim empty rows from end of each cell array
MTloc(all(cellfun(@isempty,MTloc(:,1:2)),2),:) = [];
V1loc(all(cellfun(@isempty,V1loc(:,1:2)),2),:) = [];
V1_fixloc(all(cellfun(@isempty,V1_fixloc(:,1:2)),2),:) = [];
contrast(all(cellfun(@isempty,contrast(:,1:2)),2),:) = [];
sup(all(cellfun(@isempty,sup(:,1:2)),2),:) = [];
sum(all(cellfun(@isempty,sum(:,1:2)),2),:) = [];
ftap(all(cellfun(@isempty,ftap(:,1:2)),2),:) = [];

bl_MTloc(all(cellfun(@isempty,bl_MTloc(:,1:2)),2),:) = [];
bl_V1loc(all(cellfun(@isempty,bl_V1loc(:,1:2)),2),:) = [];
bl_V1_fixloc(all(cellfun(@isempty,bl_V1_fixloc(:,1:2)),2),:) = [];
bl_contrast(all(cellfun(@isempty,bl_contrast(:,1:2)),2),:) = [];
bl_sup(all(cellfun(@isempty,bl_sup(:,1:2)),2),:) = [];
bl_sum(all(cellfun(@isempty,bl_sum(:,1:2)),2),:) = [];
bl_ftap(all(cellfun(@isempty,bl_ftap(:,1:2)),2),:) = [];


%% write cell arrays to xls file spreadsheets

success = zeros(12,1); % preallocate
if onlyGood
    xlsFilename = fullfile(save_dir,'fMRI_EyeTracking_Data_Tables_OnlyGood.xlsx'); % name file for only subjects with good data
else
    xlsFilename = fullfile(save_dir,'fMRI_EyeTracking_Data_Tables.xlsx'); % name file
end

success(1) = xlswrite(xlsFilename,MTloc,'MT Loc');
success(2) = xlswrite(xlsFilename,V1loc,'V1 Loc');
success(3) = xlswrite(xlsFilename,V1_fixloc,'V1_fix Loc');
success(4) = xlswrite(xlsFilename,contrast,'Contrast');
success(5) = xlswrite(xlsFilename,sup,'Suppression');
success(6) = xlswrite(xlsFilename,sum,'Summation');
success(7) = xlswrite(xlsFilename,ftap,'ftap');

success(8) = xlswrite(xlsFilename,bl_MTloc,'block MT Loc');
success(9) = xlswrite(xlsFilename,bl_V1loc,'block V1 Loc');
success(10) = xlswrite(xlsFilename,bl_V1_fixloc,'block V1_fix Loc');
success(11) = xlswrite(xlsFilename,bl_contrast,'block Contrast');
success(12) = xlswrite(xlsFilename,bl_sup,'block Suppression');
success(13) = xlswrite(xlsFilename,bl_sum,'block Summation');
success(14) = xlswrite(xlsFilename,bl_ftap,'block ftap');

