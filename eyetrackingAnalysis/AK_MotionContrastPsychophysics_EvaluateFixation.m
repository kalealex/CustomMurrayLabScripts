function [ subjectData ] = AK_MotionContrastPsychophysics_EvaluateFixation( subjectCodeCellStr, figuresYesNo, eventPatchesYesNo )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Check arguments
if nargin > 0 && ischar(subjectCodeCellStr)
    subjectCodeCellStr = {subjectCodeCellStr}; % make a single string input into a cell
end
if nargin < 1 || ~iscell(subjectCodeCellStr)
    error('AK_MotionContrastPsychophysics_EvaluateFixation requires subject code[s] (cell array of strings) as an argument');
end
if nargin < 3;
   eventPatchesYesNo = false;
end
if nargin < 2;
   figuresYesNo = false;
end

%% enter display characteristics for calculating visual angle, pixels per degree, and a fixation threshold

display.dist = 66; % distance to screen in cm
display.width = 36.5; % width of display in cm
display.height = 27; % height of display in cm
display.resolution = [799 599]; % number of pixels for screen width, height
display.center = display.resolution./2; % screen center in pixel coordinates
display.angles = [2*atand((display.width/2)/display.dist) 2*atand((display.height/2)/display.dist)]; % visual angle of the screen width, height
pixelsPerDegree = display.resolution./display.angles; % pixels per degree visual angle for X dimension and Y dimension

fixationRadiusThreshold = 1; % set threshold for maintaining fixation (degrees visual angle); used to be .5
fixationMinDur = 100; % fixations must be at least this many milliseconds
stimCoords = [0;0]; % coordinates of stimuli on screen (x & y by number of stim)
% trialsPerChunk_DC = 8; % number of trials to include in each chunk of drift correction (should be approximately the same as one block in fMRI to be consistent with other analysis)
minMillisecondsPerChunk = 10006; % should be at least 10 seconds per chunk of eye tracking data

%% establish directory by subjects and fMRI sessions

% set up root dir for both kinds of study
lz_top_dir = 'L:\MurrayLab\Lorazepam\Data';
gaba_top_dir = 'L:\MurrayLab\ASD\Data';
% designate folder names for invidual participants as well as lorazepam subfolders 
subj_dirs = subjectCodeCellStr; 
lz_session_dirs = {'Session 1','Session 2'};
% regular expressions for each kind of study*experiment's .mat file (these will
% be used to get the correspondence between event codes sent to the eye
% tracker and stimulus conditions)
motion_Str = '_motion_staircase'; % both lz and gaba
contrast_gabaStr = '_contrast_staircase'; % contrast detection in gaba
contrast_lzStr = '_contrast_staircase_condition'; % contrast detection/ discrimination in lz
% eyetracking name regular expressions unique to each experiment (regardless of study) experiment 
motion_eye_WildcardStr = '*m*.asc';
contrast_eye_WildcardStr = '*c*.asc';
% list of condition names corresponding to condition numbers in behavioral data files
motionCondList = {'3% contrast, 2 degrees','3% contrast, 12 degrees','98% contrast, 2 degrees','98% contrast, 12 degrees','3% contrast, 1 degree','98% contrast, 1 degree','catch trials'};
contrast_gabaCondList = {'vertical','horizontal','catch trials'};
contrast_lzCondList = {'staircase1','staircase2','staircase3','staircase4','catch trials'};%{'detection (0%)','discrimination (25%)','catch trials'};
contrast_lzExpList = {'detection', 'discrimination'};

% date of asc file parsing code
ascParseCode = dir(fullfile(gaba_top_dir,'Analysis_scripts','AMK_Functions','Eyetracking_Pipeline','AK_GABAinASD_ascfileParse.m'));
ascParseDate = datetime(ascParseCode.date);

%% cycle through subjects getting summary statistics per condition for each subject

% preallocate output struct
name = '';
directory = cell(1);
experiment(1:2) = struct;
subjectData(1:length(subj_dirs)) = struct('name',name,'directory',directory,'experiment',experiment);

for iSubj = 1:length(subj_dirs)
    % determine study based on subject code
    if ~isempty(regexp(subj_dirs{iSubj},regexptranslate('wildcard','*G*'),'once')) % is gaba participant?
        study = 'gaba';
    elseif  ~isempty(regexp(subj_dirs{iSubj},regexptranslate('wildcard','L*'),'once')) % is lorazepam participant?
        study = 'lz';
    else
        disp('You may need to modify this function in order to handle the directory structure of a new study');
        disp('or to deal with a new variant of subject codes within the GABA in ASD or Lorazepam studies.');
        error(['Subject code ' subj_dirs{iSubj} ' was not recognized as belonging to either the GABA in ASD or Lorazepam studies.'])
    end
    % store subject code in data structure
    subjectData(iSubj).name = subj_dirs{iSubj};
    % use study type as a switch in data analysis
    switch study
        case 'gaba'
            % construct working directory
            use_dir = fullfile(gaba_top_dir,subj_dirs{iSubj},'Psychophysics');
            % store directory
            subjectData(iSubj).directory = use_dir;
            %% deal with motion data:
            % set expIdx and store experiment name
            expIdx = 1;
            subjectData(iSubj).experiment(expIdx).name = 'motion staircase';
            % number of events written to eye tracker during full run
            eventsExpected = [410 590]; % could be one of two values (earlier version had only 5 conditions and thus fewer trials
            % set up run counter (runs for a session may span multiple .asc files if the experiment was stopped between runs)
            thisRun = 0;
            % set up mat file number counter (may be distinct from run count if there were bad runs)
            thisMatFileN = 1; % start checking at one
            lastMatFileN = 0; % dummy value
            % generate structure containing information about asc files
            asc_list = dir(fullfile(use_dir,motion_eye_WildcardStr));
            asc_list = AK_sortStruct(asc_list,1,[subj_dirs{iSubj} 'm']); % sort by name
            if ~isempty(asc_list)
                for iA = 1:length(asc_list); % cycle through asc files
                    % load data
                    clear block mat_name mat_date
                    mat_name = [fullfile(use_dir,asc_list(iA).name(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
                    if ~exist(mat_name,'file') % check whether or not function output has already been saved
                        disp(['parsing ' fullfile(use_dir,asc_list(iA).name)]) % message
                        block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(use_dir,asc_list(iA).name)); % parse asc file
                        save(mat_name,'block'); % save function output
                    else
                        mat_date = dir(mat_name);
                        mat_date = datetime(mat_date.date);
                        if mat_date <= ascParseDate
                            disp(['parsing ' fullfile(use_dir,asc_list(iA).name)]) % message
                            block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(use_dir,asc_list(iA).name)); % parse asc file
                            save(mat_name,'block'); % save function output
                        else
                            disp(['loading ' mat_name]) % message
                            load(mat_name,'block') % load file containing timepts array
                        end
                    end
                    for iB = 1:length(block); % cycle through blocks of trials
                        % only proceed if there are a full set events and if time is continuous
                        if any(length(block(iB).events(:,1)) == eventsExpected) && ~any((cellfun(@str2double,block(iB).timepts(3:end,1))-cellfun(@str2double,block(iB).timepts(2:end-1,1)))~=1)
                            % clear variables
                            clear time Xpixel Ypixel Xdva Ydva dist angle events eventTimes

                            % increase run count
                            thisRun = thisRun + 1;
                            
                            % store asc file name per run
                            subjectData(iSubj).experiment(expIdx).ascFilename{1,thisRun} = asc_list(iA).name;

                            % extract relevant data from block struct
                            time = cellfun(@str2double,block(iB).timepts(2:end,1)); % set time vector
                            Xpixel = cellfun(@str2double,block(iB).timepts(2:end,2)); % set vector of X position in pixels
                            Ypixel = cellfun(@str2double,block(iB).timepts(2:end,3)); % set vector of Y position in pixels
                            Xdva = (Xpixel-display.center(1))./pixelsPerDegree(1); % set vector of X position as degrees visual angle
                            Ydva = (Ypixel-display.center(2))./pixelsPerDegree(2); % set vector of Y position as degrees visual angle
                            dist = sqrt(Xdva.^2 + Ydva.^2); % calculate vector of distance from fixation as degrees visual angle
                            angle = atan2d(Ydva,Xdva); % calculate vector of polar angles at each timeptpSize = cellfun(@str2double,block(iB).timepts(2:end,4)); % set pupil size vector

                            events = block(iB).events(2:end,2); % store list of events as cell array of strings
                            eventTimes = cellfun(@str2double,block(iB).events(2:end,1)); % set vector of event times
                            events(cellfun(@isempty,events)) = {'no event'}; % fill in empty events

                            % find number of blinks and idx into time array where there are blinks
                            clear blink isBlink
                            blink.N = length(block(iB).blink(:,1))-1;
                            isBlink = zeros(length(time),1);
                            for bl = 2:length(block(iB).blink(:,1)); % cycle through saccades
                                % get saccade durations, amplitudes, and velocities
                                blink.Stimes(bl-1) = str2double(block(iB).blink{bl,1});
                                blink.Etimes(bl-1) = str2double(block(iB).blink{bl,2});
                                blink.Durs(bl-1) = str2double(block(iB).blink{bl,3});
                                % find an index in time for each blink
                                isBlink(blink.Stimes(bl-1)<= time & time<=blink.Etimes(bl-1)) = 1;
                            end

                            % find number of saccades and idx into time array where there are saccades
                            clear sacc isSacc
                            sacc.N = length(block(iB).sacc(:,1))-1;
                            isSacc = zeros(length(time),1);
                            if sacc.N>0 % is there at least one saccade?
                                for sa = 2:length(block(iB).sacc(:,1)); % cycle through saccades
                                    % get saccade start times, end times, durations, amplitudes, and velocities
                                    sacc.Stimes(sa-1) = str2double(block(iB).sacc{sa,1});
                                    sacc.Etimes(sa-1) = str2double(block(iB).sacc{sa,2});
                                    sacc.Durs(sa-1) = str2double(block(iB).sacc{sa,3}); 
                                    sacc.Amps(sa-1) = str2double(block(iB).sacc{sa,8});
                                    sacc.Velo(sa-1) = str2double(block(iB).sacc{sa,9});
                                    % find an index in time for each saccade
                                    isSacc(sacc.Stimes(sa-1)<= time & time<=sacc.Etimes(sa-1)) = 1;
                                end
                            else % if there are no saccades
                                sacc.Stimes = nan;
                                sacc.Etimes = nan;
                                sacc.Durs = nan; 
                                sacc.Amps = nan;
                                sacc.Velo = nan;
                            end

                            % determine index in time for points which are NOT saccades or blinks
                            clear fullTimeIndex NoSaccadeTimeIndex
                            fullTimeIndex = (1:length(time))';
                            if isfield(sacc,'timeIndex')
                                NoSaccadeTimeIndex = setdiff(fullTimeIndex,find(isSacc));
                            else
                                NoSaccadeTimeIndex = fullTimeIndex;
                            end

                            % differentiate saccades with blinks included
                            clear BSidx
                            BSidx = isBlink+isSacc;
                            if length(find(BSidx==0)) > 1
                                for iSa = 1:sacc.N
                                    if ~isnan(sacc.Stimes(iSa)) && ~isnan(sacc.Etimes(iSa))
                                        clear tempIdxS tempIdxE
                                        tempIdxS = find(time==sacc.Stimes(iSa));
                                        tempIdxE = find(time==sacc.Etimes(iSa));
                                        BSidx(tempIdxS:tempIdxE) = BSidx(tempIdxS:tempIdxE)*max(BSidx(tempIdxS:tempIdxE)); % distinguish blinks from saccades
                                    end
                                end
                                BSidx(BSidx>1) = nan; % use to filter blinks from saccade count (eyes open = 0, sacc = 1, blink = nan)
                            end

                            % load behavioral data and match event codes to conditions:
                            clear behavioralDataFilename thisRunMatFile conditionIdx conditionList uniqueConds conditionStartTimes condTimeIdx
                            % determine next mat file number, generate
                            % full .mat file name, make sure file exists
                            behavioralDataFilename = [subj_dirs{iSubj} motion_Str num2str(thisMatFileN) '.mat'];
                            thisRunMatFile = fullfile(use_dir,behavioralDataFilename);
                            while exist(thisRunMatFile,'file') == 0 || thisMatFileN == lastMatFileN
                                thisMatFileN = thisMatFileN + 1;
                                behavioralDataFilename = [subj_dirs{iSubj} motion_Str num2str(thisMatFileN) '.mat'];
                                thisRunMatFile = fullfile(use_dir,behavioralDataFilename);
                            end
                            lastMatFileN = thisMatFileN; % in preparation for next run
                            % load
                            load(thisRunMatFile);
                            % match up beginnings of trials (cue times) with condition numbers
                            conditionIdx = runData.conditionOrder; % order of numbered conditions
                            if max(conditionIdx) == 5 % for early version of ssMotion with only 5 conds
                                conditionIdx(conditionIdx == 5) = 7; % code last condition as catch trials
                            end
                            conditionList = motionCondList(conditionIdx); % list of conditions by name, rather than number
%                             uniqueConds = unique(conditionList);
                            uniqueConds = motionCondList(ismember(motionCondList,conditionList)); % preserve desired order but use only conditions which were run
                            conditionStartTimes = eventTimes(~cellfun(@isempty,regexp(events,regexptranslate('wildcard','*Trial*_cue_on'))));
                            condTimeIdx = arrayfun(@(x) find(time == x, 1),conditionStartTimes);
                            % store behavioral data filename and run name
                            subjectData(iSubj).experiment(expIdx).behavioralDataFilename{1,thisRun} = behavioralDataFilename;
                            subjectData(iSubj).experiment(expIdx).runName{1,thisRun} = ['run' num2str(thisRun)];

                            % NEW DRIFT CORRECTION: linear transformations optimized for each trial based on fixation coordinates
                            clear fixations chunkBoundIdx
                            % define what counts as a fixation
                            fixations = AK_defineFixations([Xdva,Ydva,time],fixationRadiusThreshold,fixationMinDur);

                            % separate fixations into chunks by trial and drift correct within each chunk
                            dcX = nan(size(Xdva)); % preallocate variables
                            dcY = nan(size(Ydva));
                            dcIdx = 1;
                            nodcTrial = zeros(length(conditionList),1); % logical to flag trials that receive no drift correction
                            % build chunks as sets of trials that are at least minMillisecondsPerChunk in duration
                            iTrial = 1;
                            chunkIdx = 1;
                            chunkDur = 0;
                            trialChunkIdx = nan(size(conditionStartTimes)); % preallocate as maximum size and trim later
                            while (iTrial + 1) <= length(conditionStartTimes)
                                % add up trials until duration exceeds minMillisecondsPerChunk
                                while (iTrial + 1) <= length(conditionStartTimes) && chunkDur < minMillisecondsPerChunk
                                    chunkDur = chunkDur + (conditionStartTimes(iTrial + 1) - conditionStartTimes(iTrial));
                                    iTrial = iTrial + 1;
                                end
                                % store trial index of chunk boundaries
                                trialChunkIdx(chunkIdx) = iTrial;
                                chunkIdx = chunkIdx + 1; % increment index
                                chunkDur = 0; % reset chunkDur for next chunk
                            end
                            trialChunkIdx(isnan(trialChunkIdx)) = []; % trim
                            chunkBoundIdx = [condTimeIdx(1);condTimeIdx(trialChunkIdx);length(time)]; % beginning and end indices in time of each chunk
                            for iChunk = 1:length(chunkBoundIdx)-1
                                %%%test message
                                disp(['Correcting chunk ' num2str(iChunk) ': timepts ' num2str(chunkBoundIdx(iChunk)) '-' num2str(chunkBoundIdx(iChunk+1))])
                                %%%
                                useFix = zeros(length(fixations),1); % reset index of fixations to use for this chunk
                                chunkFixTime = 0; % reset counter for number of milliseconds spent fixating during this chunk
                                % determine which fixations to use in the current chunk's drift correction
                                for iFix = 1:length(fixations)
                                    if isfield(fixations,'timeIdx') && ~isempty(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))) % compare indices for fixation to indices for chunk
                                        useFix(iFix) = 1; % use this fixation for drift correction in the current chunk
                                        chunkFixTime = chunkFixTime+length(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))); % count the amount of time spent fixating within this chunk (assumes 1000hz sample rate)
                                    end
                                end
                                useFix = logical(useFix);
                                % determine whether or not to apply drift correction
                                clear chunkFix chunkData chunkDCcoords
                                if any(useFix==1) && chunkFixTime >= 500 % only assign fixation coordinates where there are at least 500ms of fixation in a chunk
                                    % fixation coordinates for this chunk (n fixations by x & y)
                                    chunkFix = vertcat(fixations(useFix).coordinates)'; 
                                    % raw coords for this chunk    
                                    chunkData = horzcat(Xdva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)),Ydva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)))';
                                else 
                                    % flag non-drift-corrected trials
                                    if iChunk > length(trialChunkIdx) % this last chunk contains all remaining trials in the run
                                        nodcTrial(trialChunkIdx(iChunk-1):length(conditionStartTimes)) = 1;
                                    elseif iChunk == 1 % this first chunk includes all trials from the first until the next chunk boudary
                                        nodcTrial(1:trialChunkIdx(iChunk)) = 1;
                                    else % note all trials in this chunk as normal 
                                        nodcTrial(trialChunkIdx(iChunk-1):trialChunkIdx(iChunk)) = 1;
                                    end
                                    % replace uncorrected data with nans
                                    chunkFix = [];
                                    chunkData = nan(2,length(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)));
                                end

                                % apply drift correction to chunk
                                chunkDCcoords = AK_driftCorrect(chunkData,chunkFix,stimCoords);
                                % stitch together drift corrected data, removing repeated data points at interior chunkBounds
                                chunkDCcoords = chunkDCcoords';
                                if iChunk==1
                                    dcX(dcIdx:length(chunkDCcoords(:,1))) = chunkDCcoords(:,1);
                                    dcY(dcIdx:length(chunkDCcoords(:,2))) = chunkDCcoords(:,2);
                                    dcIdx = dcIdx + length(chunkDCcoords(:,1));
                                else
                                    dcX(dcIdx:dcIdx + length(chunkDCcoords(2:end,1)) - 1) = chunkDCcoords(2:end,1);
                                    dcY(dcIdx:dcIdx + length(chunkDCcoords(2:end,2)) - 1) = chunkDCcoords(2:end,2);
                                    dcIdx = dcIdx + length(chunkDCcoords(2:end,1));
                                end
                            end
                            % change to polar coordinates
                            dcDist = sqrt(dcX.^2 + dcY.^2); % calculate vector of distance from fixation as degrees visual angle
                            dcAngle = atan2d(dcY,dcX); % calculate vector of polar angles at each timept

                            if figuresYesNo == 1 % add figures?
                                % create and name figure
                                clear figname patchX patchY patchC
                                figname = strcat('Distance from Fixation Over Time for',' subject:',subj_dirs{iSubj},', file:',asc_list(iA).name); % concatenate figure name
                                figure('Name',figname,'NumberTitle','off'); % generate new figure with proper name

                                if eventPatchesYesNo == 1 % add patches?
                                    % create patches
                                    patchX = nan(4,length(condTimeIdx)); % preset size of patch input X
                                    patchY = nan(4,length(condTimeIdx)); % preset size of patch input Y
                                    patchC = nan(1,length(condTimeIdx)); % preset size of patch input for color
                                    for iCST = 1:length(conditionStartTimes); % cycle through conditionStartTimes
                                        if iCST ~= length(conditionStartTimes)
                                        patchX(:,iCST) = [conditionStartTimes(iCST) conditionStartTimes(iCST) conditionStartTimes(iCST+1) conditionStartTimes(iCST+1)]'; % set up X coordinates for patches
                                        patchY(:,iCST) = [0 max(dist) max(dist) 0]'; % set up Y coordinates for patches
                                        patchC(1,iCST) = find(strcmp(conditionList(iCST),uniqueConds)==1); % set color index for patches
                                        end
                                    end

                                    % set color properties of patches
                                    colormap(gray); % set colormap
                                    caxis([1 length(uniqueConds)]); % set number of color intervals

                                    for U = 1:length(uniqueConds) % dumb for loop to make legend work
                                        patch([-100 -100 -99 -99]'+U*iA,[-100 -99 -99 -100]'+U*iA,U)
                                    end

                                    % plot distance from fixation vs time and creat legend
                                    % for stimulus events
                                    patch(patchX,patchY,patchC); % map patches representing events onto axes
                                    legend(uniqueConds);
                                    hold on
                                end

                                plot(time,dist,'-b'); % plot distance at time points

                                % manipulate figure properties
                                axis_title = strcat(subj_dirs{iSubj},': ',asc_list(iA).name); % concatenate axis title
                                title(axis_title)
                                xlabel('Time (sampled every millisecond)');
                                ylabel('Distance (degrees visual angle)');
                                try
                                    axis([min(time(time~=0)) max(time) 0 max(dist)]);
                                catch
                                    axis([min(time(time~=0)) max(time) 0 5]);
                                end
                                set(gca,'XTick',1:50000:max(time))
                            end

                            subjectData(iSubj).experiment(expIdx).conditions{1,thisRun} = uniqueConds; % document condition names per asc file

                            % preallocate for speed:
                            % for runStats table
                            clear eventInstances fixationTime noTrackTime DistM DistSD AngleM AngleSD driftCorrectedFixationTime driftCorrectedDistM driftCorrectedDistSD driftCorrectedAngleM driftCorrectedAngleSD saccadeN saccadeTime saccadeAmpM saccadeAmpSD saccadeVeloM saccadeVeloSD
                            eventInstances = zeros(length(uniqueConds)+1,1);
                            fixationTime = zeros(length(uniqueConds)+1,1);
                            noTrackTime = zeros(length(uniqueConds)+1,1);
                            DistM = zeros(length(uniqueConds)+1,1);
                            DistSD = zeros(length(uniqueConds)+1,1);
                            AngleM = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            AngleSD = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            driftCorrectedFixationTime = zeros(length(uniqueConds)+1,1);
                            nonDriftCorrectedTime = zeros(length(uniqueConds)+1,1); % added 9/13/16
                            driftCorrectedDistM = zeros(length(uniqueConds)+1,1);
                            driftCorrectedDistSD = zeros(length(uniqueConds)+1,1);
                            driftCorrectedAngleM = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            driftCorrectedAngleSD = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            blinkN = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            saccadeN = zeros(length(uniqueConds)+1,1);
                            bigSaccadeN = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            saccadeTime = zeros(length(uniqueConds)+1,1);
                            saccadeAmpM = zeros(length(uniqueConds)+1,1);
                            saccadeAmpSD = zeros(length(uniqueConds)+1,1);
                            saccadeVeloM = zeros(length(uniqueConds)+1,1);
                            saccadeVeloSD = zeros(length(uniqueConds)+1,1);
                            % stats by trial
                            clear nanTime fixTime eventDur eventDistM eventDistSD eventAngleM eventAngleSD eventN ResFixTime ResDistM ResDistSD ResAngleM ResAngleSD ResN saccN bigSaccN saccDur saccAmpM saccAmpSD saccAmpN saccVeloM saccVeloSD saccVeloN
                            nanTime = cell(length(uniqueConds),1); 
                            fixTime = cell(length(uniqueConds),1); 
                            eventDur = cell(length(uniqueConds),1); 
                            eventDistM = cell(length(uniqueConds),1); 
                            eventDistSD = cell(length(uniqueConds),1); 
                            eventAngleM = cell(length(uniqueConds),1);  % added 7/27/16
                            eventAngleSD = cell(length(uniqueConds),1);  % added 7/27/16
                            eventN = cell(length(uniqueConds),1); 
                            ResFixTime = cell(length(uniqueConds),1);  
                            ResDistM = cell(length(uniqueConds),1); 
                            ResDistSD = cell(length(uniqueConds),1); 
                            ResAngleM = cell(length(uniqueConds),1);  % added 7/27/16
                            ResAngleSD = cell(length(uniqueConds),1);  % added 7/27/16
                            ResN = cell(length(uniqueConds),1); 
                            nodcTime = cell(length(uniqueConds),1); % added 9/13/16
                            blnkN = cell(length(uniqueConds),1);  % added 7/27/16
                            saccN = cell(length(uniqueConds),1);  
                            bigSaccN = cell(length(uniqueConds),1);  % added 7/27/16
                            saccDur = cell(length(uniqueConds),1);  
                            saccAmpM = cell(length(uniqueConds),1);  
                            saccAmpSD = cell(length(uniqueConds),1); 
                            saccAmpN = cell(length(uniqueConds),1); 
                            saccVeloM = cell(length(uniqueConds),1); 
                            saccVeloSD = cell(length(uniqueConds),1); 
                            saccVeloN = cell(length(uniqueConds),1); 

                            % extract statistics per condition
                            for iU = 1:length(uniqueConds); % cycle through unique events
                                % find times corresponding to event start and end
                                clear eventIndex startPos eventStart eventEnd isnotDC
                                eventIndex = strcmp(conditionList,uniqueConds(iU)); % create event index
                                startPos = find(eventIndex==1); % create vector of positions in eventTime when event starts
                                eventStart = conditionStartTimes(startPos); % event starts
                                eventEnd = nan(size(startPos)); % preallocate
                                if startPos(end) == length(conditionStartTimes)
                                    eventEnd(1:length(startPos)-1) = conditionStartTimes(startPos(1:length(startPos)-1)+1); % all events but last event end at beginning of next event
                                    eventEnd(length(startPos)) = time(end); % last event ends at end of time for this run
                                else
                                    eventEnd = conditionStartTimes(startPos+1); % each event ends at beginning of next event
                                end
                                isnotDC = logical(nodcTrial(startPos)); % logical index of which trials are not drift corrected

                                % preallocate for speed
                                clear pos 
                                pos.fixTime = cell(length(eventStart));
                                pos.time = cell(length(eventStart));
                                pos.nanTime = cell(length(eventStart));
                                pos.ResFixTime = cell(length(eventStart));
                                pos.saccSTime = cell(length(eventStart));
                                pos.saccETime = cell(length(eventStart));
                                pos.saccBlink = cell(length(eventStart));
                                pos.saccSmall = cell(length(eventStart));
                                pos.saccSaETime = cell(length(eventStart));

                                % find values of interest at each instance of event
                                for iI = 1:length(eventStart); % cycle through instances of event
                                    % raw data stats:
                                    eventDur{iU}(iI) = eventEnd(iI)-eventStart(iI); % calculate duration of each instance of event
                                    pos.time{iI} = find(eventStart(iI)<=time(1:end) & time(1:end)<eventEnd(iI)); % find positions in time during instance of event
                                    pos.nanTime{iI} = find(isnan(dist(pos.time{iI})) == 1); % find positions in where dist == nan during event
                                    nanTime{iU}(iI) = length(pos.nanTime{iI}); % calculate number of timepoints during which there is no data on each instance of event 
                                    pos.fixTime{iI} = find(dist(pos.time{iI}) < fixationRadiusThreshold); % find positions in fixTime during event
                                    fixTime{iU}(iI) = length(pos.fixTime{iI}); % calculate number of timepoints during which subject fixated on each instance of event
                                    eventDistM{iU}(iI) = nanmean(dist(pos.time{iI})); % find mean of distance measures of eye from fixation for timepoints during instance of event
                                    eventDistSD{iU}(iI) = nanstd(dist(pos.time{iI})); % find sd of distance measures of eye from fixation for timepoints during instance of event
                                    eventAngleM{iU}(iI) = nanmean(angle(pos.time{iI})); % find mean of angle measures of eye from fixation for timepoints during instance of event
                                    eventAngleSD{iU}(iI) = nanstd(angle(pos.time{iI})); % find sd of angle measures of eye from fixation for timepoints during instance of event
                                    eventN{iU}(iI) = length(dist(pos.time{iI})); % find number of distance measures for timepoints during instance of event
                                    % drift corrected data stats:
                                    if isnotDC(iI)
                                        pos.ResFixTime{iI} = []; % find positions in time during event where residualDist < fixThreshold
                                        ResFixTime{iU}(iI) = 0; % calculate number of timepoints during which subject fixated on each instance of event
                                        ResDistM{iU}(iI) = nan; % find mean of residuals from linear model for timepoints during instance of event
                                        ResDistSD{iU}(iI) = nan; % find sd of residuals from linear model for timepoints during instance of event
                                        ResAngleM{iU}(iI) = nan; % find mean of residuals from linear model for timepoints during instance of event
                                        ResAngleSD{iU}(iI) = nan; % find sd of residuals from linear model for timepoints during instance of event
                                        ResN{iU}(iI) = 0; % find number of residuals for timepoints during instance of event
                                        nodcTime{iU}(iI) = eventEnd(iI)-eventStart(iI); % amount of time where drift correction was not applied
                                    else % only use stats for drift corrected trials
                                        pos.ResFixTime{iI} = find(abs(dcDist(pos.time{iI})) < fixationRadiusThreshold); % find positions in time during event where residualDist < fixThreshold
                                        ResFixTime{iU}(iI) = length(pos.ResFixTime{iI}); % calculate number of timepoints during which subject fixated on each instance of event
                                        ResDistM{iU}(iI) = nanmean(dcDist(pos.time{iI})); % find mean of residuals from linear model for timepoints during instance of event
                                        ResDistSD{iU}(iI) = nanstd(dcDist(pos.time{iI})); % find sd of residuals from linear model for timepoints during instance of event
                                        ResAngleM{iU}(iI) = nanmean(dcAngle(pos.time{iI})); % find mean of residuals from linear model for timepoints during instance of event
                                        ResAngleSD{iU}(iI) = nanstd(dcAngle(pos.time{iI})); % find sd of residuals from linear model for timepoints during instance of event
                                        ResN{iU}(iI) = length(dcDist(pos.time{iI})); % find number of residuals for timepoints during instance of event
                                        nodcTime{iU}(iI) = 0; % amount of time where drift correction was not applied
                                    end
                                    % saccade and blink stats:
                                    pos.saccSTime{iI} = find(eventStart(iI)<=sacc.Stimes(1:end) & sacc.Stimes(1:end)<eventEnd(iI)); % find positions in sacc.Stimes during event
                                    pos.saccETime{iI} = find(eventStart(iI)<=sacc.Etimes(1:end) & sacc.Etimes(1:end)<eventEnd(iI)); % find positions in sacc.Etimes during event
                                    % add blink count
                                    pos.saccBlink{iI} = pos.saccSTime{iI}(isnan(BSidx(ismember(time,sacc.Stimes(pos.saccSTime{iI}))))); % find positions in sacc.Stimes where there was a blink
                                    blnkN{iU}(iI) = length(pos.saccBlink{iI}); % find number of blinks during instance of event
                                    % remove blinks from saccade count
                                    saccN{iU}(iI) = length(setdiff(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes(pos.saccBlink{iI}))); % find number of saccades (that are not blinks) during instance of event
                                    % create separate saccade count that excludes saccades w/ amp>fixThreshold
                                    pos.saccSmall{iI} = pos.saccSTime{iI}(ismember(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes(sacc.Amps < fixationRadiusThreshold))==1); % find positions in sacc.Stimes where the saccade amplitude < fixThreshold
                                    bigSaccN{iU}(iI) = length(setdiff(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes([pos.saccBlink{iI} pos.saccSmall{iI}])));
                                    pos.saccSaETime{iI} = intersect(pos.saccSTime{iI},pos.saccETime{iI}); % find positions in array of saccades where saccades start and end within instance of event
                                    saccDur{iU}(iI) = sum(sacc.Durs(pos.saccSaETime{iI})); % find the sum of saccade durations which occur within instance of event
                                    if length(pos.saccSTime{iI}) >= length(pos.saccETime{iI}) % find longer pos.sacc*Time array
                                        for st = 1:length(pos.saccSTime{iI}); %cycle through saccades starting during instance of event
                                            % add to saccDur(iI) the within-instance duration of
                                            % saccades starting, but not ending, in event instance 
                                            if st <= length(pos.saccSTime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccSTime{iI}(st)))
                                                saccDur{iU}(iI) = saccDur{iU}(iI) + (eventEnd(iI) - sacc.Stimes(pos.saccSTime{iI}(st)));
                                            % add to saccDur(iI) the within-instance duration of
                                            % saccades ending, but not starting, in event instance
                                            elseif st <= length(pos.saccETime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccETime{iI}(st)))
                                                saccDur{iU}(iI) = saccDur{iU}(iI) + (sacc.Etimes(pos.saccETime{iI}(st)) - eventStart(iI));
                                            end
                                        end
                                    elseif length(pos.saccSTime{iI}) < length(pos.saccETime{iI}) % find longer pos.sacc*Time array
                                        for st = 1:length(pos.saccETime{iI}); %cycle through saccades ending during instance of event
                                            % add to saccDur(i) the within-instance duration of
                                            % saccades starting, but not ending, in event instance 
                                            if st <= length(pos.saccSTime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccSTime{iI}(st)))
                                                saccDur{iU}(iI) = saccDur{iU}(iI) + (eventEnd(iI) - sacc.Stimes(pos.saccSTime{iI}(st)));
                                            % add to saccDur(i) the within-instance duration of
                                            % saccades ending, but not starting, in event instance
                                            elseif st <= length(pos.saccETime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccETime{iI}(st)))
                                                saccDur{iU}(iI) = saccDur{iU}(iI) + (sacc.Etimes(pos.saccETime{iI}(st)) - eventStart(iI));
                                            end
                                        end
                                    end
                                    saccAmpM{iU}(iI) = nanmean(sacc.Amps(pos.saccSTime{iI})); % find mean of amp for timepoints during instance of event
                                    saccAmpSD{iU}(iI) = nanstd(sacc.Amps(pos.saccSTime{iI})); % find sd of amp for timepoints during instance of event
                                    saccAmpN{iU}(iI) = length(sacc.Amps(pos.saccSTime{iI})); % find number of amp measures for timepoints during instance of event
                                    saccVeloM{iU}(iI) = nanmean(sacc.Velo(pos.saccSTime{iI})); % find mean of velocity for timepoints during instance of event
                                    saccVeloSD{iU}(iI) = nanstd(sacc.Velo(pos.saccSTime{iI})); % find sd of velocity for timepoints during instance of event
                                    saccVeloN{iU}(iI) = length(sacc.Velo(pos.saccSTime{iI})); % find number of velocity measures for timepoints during instance of event
                                end

                                % calculate statistics to be saved in data structure
                                if ~isempty(uniqueConds(iU)) % only assign meaningful values to data structure for actual events
                                    % document number of instances of event
                                    eventInstances(iU) = length(eventStart);

                                    % runStats:
                                    % find nan time (%) for event
                                    clear TotalNanTime TotalEventDur
                                    TotalNanTime = nansum(nanTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                    TotalEventDur = nansum(eventDur{iU}(1:end)); % total duration of all instances of event
                                    noTrackTime(iU) = TotalNanTime/TotalEventDur; % document proportion of time fixating across instances of event

                                    % find fixation time (%) for event
                                    clear TotalFixTime
                                    TotalFixTime = nansum(fixTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                    fixationTime(iU) = TotalFixTime/TotalEventDur; % document proportion of time fixating across instances of event

                                    % find distance and angle from fixation stats during event
                                    clear DistGM DistGSD AngleGM AngleGSD
                                    [DistGM, DistGSD] = AK_grandSD(eventDistM{iU},eventDistSD{iU},eventN{iU}); % find mean and SD for all instances of this event
                                    DistM(iU) = DistGM; % document function output
                                    DistSD(iU) = DistGSD; % document function output
                                    [AngleGM, AngleGSD] = AK_grandSD(eventAngleM{iU},eventAngleSD{iU},eventN{iU}); % find mean and SD for all instances of this event
                                    AngleM(iU) = AngleGM; % document function output
                                    AngleSD(iU) = AngleGSD; % document function output

                                    % find drift-corrected fixation time (%) for event
                                    clear TotalResFixTime
                                    TotalResFixTime = nansum(ResFixTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                    driftCorrectedFixationTime(iU) = TotalResFixTime/TotalEventDur; % document proportion of time fixating across instances of event

                                    % find proportion of time for which drift correction was not applied
                                    clear TotalNoDCTime
                                    TotalNoDCTime = nansum(nodcTime{iU}(1:end)); % number of timepoints not drift corrected across all instances of event
                                    nonDriftCorrectedTime(iU) = TotalNoDCTime/TotalEventDur; % document proportion of time fixating across instances of event

                                    % find drift-corrected distance and angle from fixation stats during event
                                    clear ResDistGM ResDistGSD ResAngleGM ResAngleGSD
                                    [ResDistGM, ResDistGSD] = AK_grandSD(ResDistM{iU},ResDistSD{iU},ResN{iU}); % find mean and SD for all instances of this event
                                    driftCorrectedDistM(iU) = ResDistGM; % document function output
                                    driftCorrectedDistSD(iU) = ResDistGSD; % document function output
                                    [ResAngleGM, ResAngleGSD] = AK_grandSD(ResAngleM{iU},ResAngleSD{iU},ResN{iU}); % find mean and SD for all instances of this event
                                    driftCorrectedAngleM(iU) = ResAngleGM; % document function output
                                    driftCorrectedAngleSD(iU) = ResAngleGSD; % document function output

                                    % find saccade stats during event
                                    clear TotalSaccTime
                                    blinkN(iU) = nansum(blnkN{iU}(1:end)); % document sum of number of blinks across instances of event
                                    saccadeN(iU) = nansum(saccN{iU}(1:end)); % document sum of number of sacccades across instances of event
                                    bigSaccadeN(iU) = nansum(bigSaccN{iU}(1:end)); % document sum of number of sacccades outside of fixation threshold across instances of event
                                    TotalSaccTime = nansum(saccDur{iU}(1:end)); % sum saccade durations across instances of event
                                    saccadeTime(iU) = TotalSaccTime/TotalEventDur; % document proportion of time saccading across instances of event
                                    % saccade amplitude stats
                                    clear saccAmpGM saccAmpGSD
                                    [saccAmpGM, saccAmpGSD] = AK_grandSD(saccAmpM{iU},saccAmpSD{iU},saccAmpN{iU}); % find mean and SD for all instances of this event
                                    saccadeAmpM(iU) = saccAmpGM; % document function output
                                    saccadeAmpSD(iU) = saccAmpGSD; % document function output
                                    % saccade velocity stats
                                    clear saccVeloGM saccVeloGSD
                                    [saccVeloGM, saccVeloGSD] = AK_grandSD(saccVeloM{iU},saccVeloSD{iU},saccVeloN{iU}); % find mean and SD for all instances of this event
                                    saccadeVeloM(iU) = saccVeloGM; % document function output
                                    saccadeVeloSD(iU) = saccVeloGSD; % document function output

                                elseif isempty(uniqueConds(iU)) % where uniqueConds(iU) does not exist for file(iA), set stats values = []
                                    eventInstances(iU) = nan;
                                    fixationTime(iU) = nan;
                                    noTrackTime(iU) = nan;
                                    DistM(iU) = nan;
                                    DistSD(iU) = nan;
                                    AngleM(iU) = nan;
                                    AngleSD(iU) = nan;
                                    driftCorrectedFixationTime(iU) = nan;
                                    driftCorrectedDistM(iU) = nan;
                                    driftCorrectedDistSD(iU) = nan;
                                    driftCorrectedAngleM(iU) = nan;
                                    driftCorrectedAngleSD(iU) = nan;
                                    blinkN(iU) = nan;
                                    saccadeN(iU) = nan;
                                    bigSaccadeN(iU) = nan;
                                    saccadeTime(iU) = nan;
                                    saccadeAmpM(iU) = nan;
                                    saccadeAmpSD(iU) = nan;
                                    saccadeVeloM(iU) = nan;
                                    saccadeVeloSD(iU) = nan;
                                end
                            end

                            % calculate statistics across all conditions
                            eventInstances(iU+1) = []; % N/A
                            pos.ALLfixTime = find(dist < fixationRadiusThreshold); % find positions where distance < fixThreshold
                            fixationTime(iU+1) = length(pos.ALLfixTime)/length(time); % calculate fixation time (%) for whole file
                            pos.ALLnanTime = find(isnan(dist) == 1); % find positions in where dist == nan 
                            noTrackTime(iU+1) = length(pos.ALLnanTime)/length(time); % calculate time (%) where there is no tracking
                            DistM(iU+1) = nanmean(dist); % calculate mean distance from fixation
                            DistSD(iU+1) = nanstd(dist); % calculate SD of distance from fixation
                            AngleM(iU+1) = nanmean(angle); % calculate mean angle from fixation
                            AngleSD(iU+1) = nanstd(angle); % calculate SD of angle from fixation
                            pos.ALLResFixTime = find(dcDist < fixationRadiusThreshold); % find positions where drift corrected distance from fixation < fixThreshold
                            driftCorrectedFixationTime(iU+1) = length(pos.ALLResFixTime)/length(time); % calculate drift corrected fixation time (%) for whole file
                            nonDriftCorrectedTime(iU+1) = nansum(nonDriftCorrectedTime(1:iU))/length(time);
                            driftCorrectedDistM(iU+1) = nanmean(dcDist); % calculate mean drift corrected distance from fixation
                            driftCorrectedDistSD(iU+1) = nanstd(dcDist); % calculate SD of drift corrected distance from fixation
                            driftCorrectedAngleM(iU+1) = nanmean(dcAngle); % calculate mean drift corrected distance from fixation
                            driftCorrectedAngleSD(iU+1) = nanstd(dcAngle); % calculate SD of drift corrected distance from fixation
                            blinkN(iU+1) = nansum(blinkN(1:end-1)); % document total number of blinks in whole file
                            saccadeN(iU+1) =  nansum(saccadeN(1:end-1)); % document total number of saccades in whole file
                            bigSaccadeN(iU+1) =  nansum(bigSaccadeN(1:end-1)); % document total number of saccades outside of fixation threshold in whole file
                            saccadeTime(iU+1) = nansum(sacc.Durs)/length(time); % calculate time (%) where subject is in a saccade
                            saccadeAmpM(iU+1) = nanmean(sacc.Amps); % calculate the mean saccade amplitude for the whole file
                            saccadeAmpSD(iU+1) = nanstd(sacc.Amps); % calculate the SD of saccade amplitudes for the whole file
                            saccadeVeloM(iU+1) = nanmean(sacc.Velo); % calculate the mean saccade velocity for the whole file
                            saccadeVeloSD(iU+1) = nanstd(sacc.Velo); % calculate the SD of saccade velocities for the whole file

                            % store stats from asc file iA in data structure
                            subjectData(iSubj).experiment(expIdx).conditionTrialCount{1,thisRun} = eventInstances; % document the number of instances of each condition

                            % designate labels for stats tables
                            clear summaryStatsTable_condLabels statsTable_statLabels
                            runStatsTable_condLabels = [uniqueConds'; {'all'}]';
                            statsTable_statLabels = {'fixationTime','driftCorrectedFixationTime','noDriftCorrectTime','noTrackTime','saccadeTime','distanceMean','distanceSD','angleMean','angleSD','driftCorrectedDistanceMean','driftCorrectedDistanceSD','driftCorrectedAngleMean','driftCorrectedAngleSD','BlinkN','SaccadeN','bigSaccadeN','SaccadeAmplitudeMean','SaccadeAmplitudeSD','SaccadeVelocityMean','SaccadeVelocitySD'};

                            % document table of summary statistics for
                            % each condition in each run
                            subjectData(iSubj).experiment(expIdx).runStats{1,thisRun} = table(fixationTime,driftCorrectedFixationTime,nonDriftCorrectedTime,noTrackTime,saccadeTime,DistM,DistSD,AngleM,AngleSD,driftCorrectedDistM,driftCorrectedDistSD,driftCorrectedAngleM,driftCorrectedAngleSD,blinkN,saccadeN,bigSaccadeN,saccadeAmpM,saccadeAmpSD,saccadeVeloM,saccadeVeloSD,'RowNames',runStatsTable_condLabels,'VariableNames',statsTable_statLabels);

                            % assess data quality for asc file
                            % data is bad quality if all data is missing,
                            % there is no data for at least 40% of the run,
                            % or the standard deviation of the drift
                            % corrected distrance from fixation outside of
                            % saccade events is greater than 1.5 degrees
                            % visual angle
                            if mean(isnan(dist)) == 1 || noTrackTime(end) > .40 || nanstd(dcDist(NoSaccadeTimeIndex)) > 1.5 
                                subjectData(iSubj).experiment(expIdx).quality{1,thisRun} = 0;
                            else
                                subjectData(iSubj).experiment(expIdx).quality{1,thisRun} = 1;
                            end                         
                        end
                    end
                end
            else
                % fill in fields with nans because there are no asc
                % files for this session directory
                subjectData(iSubj).experiment(expIdx).ascFilename = nan; 
                subjectData(iSubj).experiment(expIdx).behavioralDataFilename = nan; 
                subjectData(iSubj).experiment(expIdx).conditions = nan;
                subjectData(iSubj).experiment(expIdx).conditionTrialCount = nan;
                subjectData(iSubj).experiment(expIdx).runStats = nan;
                subjectData(iSubj).experiment(expIdx).quality = 'could not find asc files';
                disp(['No ' motion_eye_WildcardStr ' files for subject ' subj_dirs{iSubj}]);
            end
            %% deal with contrast data:
            % set expIdx and store experiment name
            expIdx = 2;
            subjectData(iSubj).experiment(expIdx).name = 'contrast staircase';
            % number of events written to eye tracker during full run
            eventsExpected = [236 248]; % odd numbered runs have a different number of trials from even numbered runs
            % set up run counter (runs for a session may span multiple .asc files if the experiment was stopped between runs)
            thisRun = 0;
            % generate structure containing information about asc files
            asc_list = dir(fullfile(use_dir,contrast_eye_WildcardStr)); 
            asc_list = AK_sortStruct(asc_list,1,[subj_dirs{iSubj} 'c']); % sort by name
            if ~isempty(asc_list)
                for iA = 1:length(asc_list); % cycle through asc files
                    % load data
                    clear block mat_name mat_date
                    mat_name = [fullfile(use_dir,asc_list(iA).name(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
                    if ~exist(mat_name,'file') % check whether or not function output has already been saved
                        disp(['parsing ' fullfile(use_dir,asc_list(iA).name)]) % message
                        block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(use_dir,asc_list(iA).name)); % parse asc file
                        save(mat_name,'block'); % save function output
                    else
                        mat_date = dir(mat_name);
                        mat_date = datetime(mat_date.date);
                        if mat_date <= ascParseDate
                            disp(['parsing ' fullfile(use_dir,asc_list(iA).name)]) % message
                            block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(use_dir,asc_list(iA).name)); % parse asc file
                            save(mat_name,'block'); % save function output
                        else
                            disp(['loading ' mat_name]) % message
                            load(mat_name,'block') % load file containing timepts array
                        end
                    end
                    % determine number of events in each block
                    clear blockEventCount
                    blockEventCount = AK_structfun(@length,block,'events','(:,1)');
                    for iB = 1:length(block); % cycle through blocks of trials
                        % only proceed if there are a full set events and
                        % if time is continuous and if both runs within the
                        % staircase ran to completion (to catch bug for G106, where experiment was terminated mid-staircase)
                        if ~any((cellfun(@str2double,block(iB).timepts(3:end,1))-cellfun(@str2double,block(iB).timepts(2:end-1,1)))~=1) && ...
                                ((blockEventCount(iB) == eventsExpected(1) && any(blockEventCount(iB:end) == eventsExpected(2))) || blockEventCount(iB) == eventsExpected(2))
                            % clear variables
                            clear time Xpixel Ypixel Xdva Ydva dist angle events eventTimes

                            % increase run count
                            thisRun = thisRun + 1;
                            % deal with bug for G114,
                            % where run1 edf did not write fully
                            if blockEventCount(iB) == eventsExpected(2) && mod(thisRun,2) == 1 % if the code is about to process an even numbered run as an odd numbered run 
                                % fill in first run output with nans
                                subjectData(iSubj).experiment(expIdx).ascFilename{1,thisRun} = nan;
                                subjectData(iSubj).experiment(expIdx).behavioralDataFilename{1,thisRun} = nan;
                                subjectData(iSubj).experiment(expIdx).runName{1,thisRun} = nan;
                                subjectData(iSubj).experiment(expIdx).conditions{1,thisRun} = nan;
                                subjectData(iSubj).experiment(expIdx).conditionTrialCount{1,thisRun} = nan;
                                subjectData(iSubj).experiment(expIdx).runStats{1,thisRun} = nan;
                                subjectData(iSubj).experiment(expIdx).quality{1,thisRun} = 'eyetracking data in asc file was incomplete';
                                % increment run count
                                thisRun = thisRun + 1;
                            end
                            % store asc file name per run
                            subjectData(iSubj).experiment(expIdx).ascFilename{1,thisRun} = asc_list(iA).name;

                            % extract relevant data from block struct
                            time = cellfun(@str2double,block(iB).timepts(2:end,1)); % set time vector
                            Xpixel = cellfun(@str2double,block(iB).timepts(2:end,2)); % set vector of X position in pixels
                            Ypixel = cellfun(@str2double,block(iB).timepts(2:end,3)); % set vector of Y position in pixels
                            Xdva = (Xpixel-display.center(1))./pixelsPerDegree(1); % set vector of X position as degrees visual angle
                            Ydva = (Ypixel-display.center(2))./pixelsPerDegree(2); % set vector of Y position as degrees visual angle
                            dist = sqrt(Xdva.^2 + Ydva.^2); % calculate vector of distance from fixation as degrees visual angle
                            angle = atan2d(Ydva,Xdva); % calculate vector of polar angles at each timeptpSize = cellfun(@str2double,block(iB).timepts(2:end,4)); % set pupil size vector

                            events = block(iB).events(2:end,2); % store list of events as cell array of strings
                            eventTimes = cellfun(@str2double,block(iB).events(2:end,1)); % set vector of event times
                            events(cellfun(@isempty,events)) = {'no event'}; % fill in empty events

                            % find number of blinks and idx into time array where there are blinks
                            clear blink isBlink
                            blink.N = length(block(iB).blink(:,1))-1;
                            isBlink = zeros(length(time),1);
                            for bl = 2:length(block(iB).blink(:,1)); % cycle through saccades
                                % get saccade durations, amplitudes, and velocities
                                blink.Stimes(bl-1) = str2double(block(iB).blink{bl,1});
                                blink.Etimes(bl-1) = str2double(block(iB).blink{bl,2});
                                blink.Durs(bl-1) = str2double(block(iB).blink{bl,3});
                                % find an index in time for each blink
                                isBlink(blink.Stimes(bl-1)<= time & time<=blink.Etimes(bl-1)) = 1;
                            end

                            % find number of saccades and idx into time array where there are saccades
                            clear sacc isSacc
                            sacc.N = length(block(iB).sacc(:,1))-1;
                            isSacc = zeros(length(time),1);
                            if sacc.N>0 % is there at least one saccade?
                                for sa = 2:length(block(iB).sacc(:,1)); % cycle through saccades
                                    % get saccade start times, end times, durations, amplitudes, and velocities
                                    sacc.Stimes(sa-1) = str2double(block(iB).sacc{sa,1});
                                    sacc.Etimes(sa-1) = str2double(block(iB).sacc{sa,2});
                                    sacc.Durs(sa-1) = str2double(block(iB).sacc{sa,3}); 
                                    sacc.Amps(sa-1) = str2double(block(iB).sacc{sa,8});
                                    sacc.Velo(sa-1) = str2double(block(iB).sacc{sa,9});
                                    % find an index in time for each saccade
                                    isSacc(sacc.Stimes(sa-1)<= time & time<=sacc.Etimes(sa-1)) = 1;
                                end
                            else % if there are no saccades
                                sacc.Stimes = nan;
                                sacc.Etimes = nan;
                                sacc.Durs = nan; 
                                sacc.Amps = nan;
                                sacc.Velo = nan;
                            end

                            % determine index in time for points which are NOT saccades or blinks
                            clear fullTimeIndex NoSaccadeTimeIndex
                            fullTimeIndex = (1:length(time))';
                            if isfield(sacc,'timeIndex')
                                NoSaccadeTimeIndex = setdiff(fullTimeIndex,find(isSacc));
                            else
                                NoSaccadeTimeIndex = fullTimeIndex;
                            end

                            % differentiate saccades with blinks included
                            clear BSidx
                            BSidx = isBlink+isSacc;
                            if length(find(BSidx==0)) > 1
                                for iSa = 1:sacc.N
                                    if ~isnan(sacc.Stimes(iSa)) && ~isnan(sacc.Etimes(iSa))
                                        clear tempIdxS tempIdxE
                                        tempIdxS = find(time==sacc.Stimes(iSa));
                                        tempIdxE = find(time==sacc.Etimes(iSa));
                                        BSidx(tempIdxS:tempIdxE) = BSidx(tempIdxS:tempIdxE)*max(BSidx(tempIdxS:tempIdxE)); % distinguish blinks from saccades
                                    end
                                end
                                BSidx(BSidx>1) = nan; % use to filter blinks from saccade count (eyes open = 0, sacc = 1, blink = nan)
                            end

                            % load behavioral data and match event codes to conditions:
                            % generate full .mat file name and load
                            clear contrastMatfileN behavioralDataFilename thisRunMatFile conditionIdx conditionList uniqueConds conditionStartTimes condTimeIdx
                            if thisRun > 2
                                contrastMatfileN = 2;
                            else % if run 1 or 2
                                contrastMatfileN = 1;
                            end
                            behavioralDataFilename = [subj_dirs{iSubj} contrast_gabaStr num2str(contrastMatfileN) '.mat'];
                            thisRunMatFile = fullfile(use_dir,behavioralDataFilename);
                            load(thisRunMatFile);
                            % match up beginnings of trials (cue times) with condition numbers:
                            % contrast detection experiment is organized
                            % into two staircases, each with a set of two
                            % runs: the first run in each staircase has 39
                            % trials; the second run in each staircase has
                            % 41 trials
                            if thisRun == 1 || thisRun == 3
                                conditionIdx = runData.conditionOrder(1:39); % order of numbered conditions
                            elseif thisRun == 2 || thisRun == 4
                                conditionIdx = runData.conditionOrder(40:end); % order of numbered conditions
                            end
                            conditionList = contrast_gabaCondList(conditionIdx); % list of conditions by name, rather than number
%                             uniqueConds = unique(conditionList);
                            uniqueConds = contrast_gabaCondList(ismember(contrast_gabaCondList,conditionList)); % preserve desired order but use only conditions which were run
                            conditionStartTimes = eventTimes(~cellfun(@isempty,regexp(events,regexptranslate('wildcard','*Trial*_cue1_on'))));
                            condTimeIdx = arrayfun(@(x) find(time == x, 1),conditionStartTimes);
                            % store behavioral data filename and run name
                            subjectData(iSubj).experiment(expIdx).behavioralDataFilename{1,thisRun} = behavioralDataFilename;
                            subjectData(iSubj).experiment(expIdx).runName{1,thisRun} = contrast_lzExpList{contrastMatfileN};

                            % NEW DRIFT CORRECTION: linear transformations optimized for each trial based on fixation coordinates
                            clear fixations chunkBoundIdx
                            % define what counts as a fixation
                            fixations = AK_defineFixations([Xdva,Ydva,time],fixationRadiusThreshold,fixationMinDur);

                            % separate fixations into chunks by trial and drift correct within each chunk
                            dcX = nan(size(Xdva)); % preallocate variables
                            dcY = nan(size(Ydva));
                            dcIdx = 1;
                            nodcTrial = zeros(length(conditionList),1); % logical to flag trials that receive no drift correction
                            % build chunks as sets of trials that are at least minMillisecondsPerChunk in duration
                            iTrial = 1;
                            chunkIdx = 1;
                            chunkDur = 0;
                            trialChunkIdx = nan(size(conditionStartTimes)); % preallocate as maximum size and trim later
                            while (iTrial + 1) <= length(conditionStartTimes)
                                % add up trials until duration exceeds minMillisecondsPerChunk
                                while (iTrial + 1) <= length(conditionStartTimes) && chunkDur < minMillisecondsPerChunk
                                    chunkDur = chunkDur + (conditionStartTimes(iTrial + 1) - conditionStartTimes(iTrial));
                                    iTrial = iTrial + 1;
                                end
                                % store trial index of chunk boundaries
                                trialChunkIdx(chunkIdx) = iTrial;
                                chunkIdx = chunkIdx + 1; % increment index
                                chunkDur = 0; % reset chunkDur for next chunk
                            end
                            trialChunkIdx(isnan(trialChunkIdx)) = []; % trim
                            chunkBoundIdx = [condTimeIdx(1);condTimeIdx(trialChunkIdx);length(time)]; % beginning and end indices in time of each chunk
                            for iChunk = 1:length(chunkBoundIdx)-1
                                %%%test message
                                disp(['Correcting chunk ' num2str(iChunk) ': timepts ' num2str(chunkBoundIdx(iChunk)) '-' num2str(chunkBoundIdx(iChunk+1))])
                                %%%
                                useFix = zeros(length(fixations),1); % reset index of fixations to use for this chunk
                                chunkFixTime = 0; % reset counter for number of milliseconds spent fixating during this chunk
                                % determine which fixations to use in the current chunk's drift correction
                                for iFix = 1:length(fixations)
                                    if isfield(fixations,'timeIdx') && ~isempty(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))) % compare indices for fixation to indices for chunk
                                        useFix(iFix) = 1; % use this fixation for drift correction in the current chunk
                                        chunkFixTime = chunkFixTime+length(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))); % count the amount of time spent fixating within this chunk (assumes 1000hz sample rate)
                                    end
                                end
                                useFix = logical(useFix);
                                % determine whether or not to apply drift correction
                                clear chunkFix chunkData chunkDCcoords
                                if any(useFix==1) && chunkFixTime >= 500 % only assign fixation coordinates where there are at least 500ms of fixation in a chunk
                                    % fixation coordinates for this chunk (n fixations by x & y)
                                    chunkFix = vertcat(fixations(useFix).coordinates)'; 
                                    % raw coords for this chunk    
                                    chunkData = horzcat(Xdva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)),Ydva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)))';
                                else 
                                    % flag non-drift-corrected trials
                                    if iChunk > length(trialChunkIdx) % this last chunk contains all remaining trials in the run
                                        nodcTrial(trialChunkIdx(iChunk-1):length(conditionStartTimes)) = 1;
                                    elseif iChunk == 1 % this first chunk includes all trials from the first until the next chunk boudary
                                        nodcTrial(1:trialChunkIdx(iChunk)) = 1;
                                    else % note all trials in this chunk as normal 
                                        nodcTrial(trialChunkIdx(iChunk-1):trialChunkIdx(iChunk)) = 1;
                                    end
                                    % replace uncorrected data with nans
                                    chunkFix = [];
                                    chunkData = nan(2,length(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)));
                                end

                                % apply drift correction to chunk
                                chunkDCcoords = AK_driftCorrect(chunkData,chunkFix,stimCoords);
                                % stitch together drift corrected data, removing repeated data points at interior chunkBounds
                                chunkDCcoords = chunkDCcoords';
                                if iChunk==1
                                    dcX(dcIdx:length(chunkDCcoords(:,1))) = chunkDCcoords(:,1);
                                    dcY(dcIdx:length(chunkDCcoords(:,2))) = chunkDCcoords(:,2);
                                    dcIdx = dcIdx + length(chunkDCcoords(:,1));
                                else
                                    dcX(dcIdx:dcIdx + length(chunkDCcoords(2:end,1)) - 1) = chunkDCcoords(2:end,1);
                                    dcY(dcIdx:dcIdx + length(chunkDCcoords(2:end,2)) - 1) = chunkDCcoords(2:end,2);
                                    dcIdx = dcIdx + length(chunkDCcoords(2:end,1));
                                end
                            end
                            % change to polar coordinates
                            dcDist = sqrt(dcX.^2 + dcY.^2); % calculate vector of distance from fixation as degrees visual angle
                            dcAngle = atan2d(dcY,dcX); % calculate vector of polar angles at each timept

                            if figuresYesNo == 1 % add figures?
                                % create and name figure
                                clear figname patchX patchY patchC
                                figname = strcat('Distance from Fixation Over Time for',' subject:',subj_dirs{iSubj},', file:',asc_list(iA).name); % concatenate figure name
                                figure('Name',figname,'NumberTitle','off'); % generate new figure with proper name

                                if eventPatchesYesNo == 1 % add patches?
                                    % create patches
                                    patchX = nan(4,length(condTimeIdx)); % preset size of patch input X
                                    patchY = nan(4,length(condTimeIdx)); % preset size of patch input Y
                                    patchC = nan(1,length(condTimeIdx)); % preset size of patch input for color
                                    for iCST = 1:length(conditionStartTimes); % cycle through conditionStartTimes
                                        if iCST ~= length(conditionStartTimes)
                                        patchX(:,iCST) = [conditionStartTimes(iCST) conditionStartTimes(iCST) conditionStartTimes(iCST+1) conditionStartTimes(iCST+1)]'; % set up X coordinates for patches
                                        patchY(:,iCST) = [0 max(dist) max(dist) 0]'; % set up Y coordinates for patches
                                        patchC(1,iCST) = find(strcmp(conditionList(iCST),uniqueConds)==1); % set color index for patches
                                        end
                                    end

                                    % set color properties of patches
                                    colormap(gray); % set colormap
                                    caxis([1 length(uniqueConds)]); % set number of color intervals

                                    for U = 1:length(uniqueConds) % dumb for loop to make legend work
                                        patch([-100 -100 -99 -99]'+U*iA,[-100 -99 -99 -100]'+U*iA,U)
                                    end

                                    % plot distance from fixation vs time and creat legend
                                    % for stimulus events
                                    patch(patchX,patchY,patchC); % map patches representing events onto axes
                                    legend(uniqueConds);
                                    hold on
                                end

                                plot(time,dist,'-b'); % plot distance at time points

                                % manipulate figure properties
                                axis_title = strcat(subj_dirs{iSubj},': ',asc_list(iA).name); % concatenate axis title
                                title(axis_title)
                                xlabel('Time (sampled every millisecond)');
                                ylabel('Distance (degrees visual angle)');
                                try
                                    axis([min(time(time~=0)) max(time) 0 max(dist)]);
                                catch
                                    axis([min(time(time~=0)) max(time) 0 5]);
                                end
                                set(gca,'XTick',1:50000:max(time))
                            end

                            subjectData(iSubj).experiment(expIdx).conditions{1,thisRun} = uniqueConds; % document condition names per asc file

                            % preallocate for speed:
                            % for runStats table
                            clear eventInstances fixationTime noTrackTime DistM DistSD AngleM AngleSD driftCorrectedFixationTime driftCorrectedDistM driftCorrectedDistSD driftCorrectedAngleM driftCorrectedAngleSD saccadeN saccadeTime saccadeAmpM saccadeAmpSD saccadeVeloM saccadeVeloSD
                            eventInstances = zeros(length(uniqueConds)+1,1);
                            fixationTime = zeros(length(uniqueConds)+1,1);
                            noTrackTime = zeros(length(uniqueConds)+1,1);
                            DistM = zeros(length(uniqueConds)+1,1);
                            DistSD = zeros(length(uniqueConds)+1,1);
                            AngleM = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            AngleSD = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            driftCorrectedFixationTime = zeros(length(uniqueConds)+1,1);
                            nonDriftCorrectedTime = zeros(length(uniqueConds)+1,1); % added 9/13/16
                            driftCorrectedDistM = zeros(length(uniqueConds)+1,1);
                            driftCorrectedDistSD = zeros(length(uniqueConds)+1,1);
                            driftCorrectedAngleM = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            driftCorrectedAngleSD = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            blinkN = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            saccadeN = zeros(length(uniqueConds)+1,1);
                            bigSaccadeN = zeros(length(uniqueConds)+1,1); % added 7/27/16
                            saccadeTime = zeros(length(uniqueConds)+1,1);
                            saccadeAmpM = zeros(length(uniqueConds)+1,1);
                            saccadeAmpSD = zeros(length(uniqueConds)+1,1);
                            saccadeVeloM = zeros(length(uniqueConds)+1,1);
                            saccadeVeloSD = zeros(length(uniqueConds)+1,1);
                            % stats by trial
                            clear nanTime fixTime eventDur eventDistM eventDistSD eventAngleM eventAngleSD eventN ResFixTime ResDistM ResDistSD ResAngleM ResAngleSD ResN saccN bigSaccN saccDur saccAmpM saccAmpSD saccAmpN saccVeloM saccVeloSD saccVeloN
                            nanTime = cell(length(uniqueConds),1); 
                            fixTime = cell(length(uniqueConds),1); 
                            eventDur = cell(length(uniqueConds),1); 
                            eventDistM = cell(length(uniqueConds),1); 
                            eventDistSD = cell(length(uniqueConds),1); 
                            eventAngleM = cell(length(uniqueConds),1);  % added 7/27/16
                            eventAngleSD = cell(length(uniqueConds),1);  % added 7/27/16
                            eventN = cell(length(uniqueConds),1); 
                            ResFixTime = cell(length(uniqueConds),1);  
                            ResDistM = cell(length(uniqueConds),1); 
                            ResDistSD = cell(length(uniqueConds),1); 
                            ResAngleM = cell(length(uniqueConds),1);  % added 7/27/16
                            ResAngleSD = cell(length(uniqueConds),1);  % added 7/27/16
                            ResN = cell(length(uniqueConds),1); 
                            nodcTime = cell(length(uniqueConds),1); % added 9/13/16
                            blnkN = cell(length(uniqueConds),1);  % added 7/27/16
                            saccN = cell(length(uniqueConds),1);  
                            bigSaccN = cell(length(uniqueConds),1);  % added 7/27/16
                            saccDur = cell(length(uniqueConds),1);  
                            saccAmpM = cell(length(uniqueConds),1);  
                            saccAmpSD = cell(length(uniqueConds),1); 
                            saccAmpN = cell(length(uniqueConds),1); 
                            saccVeloM = cell(length(uniqueConds),1); 
                            saccVeloSD = cell(length(uniqueConds),1); 
                            saccVeloN = cell(length(uniqueConds),1); 

                            % extract statistics per condition
                            for iU = 1:length(uniqueConds); % cycle through unique events
                                % find times corresponding to event start and end
                                clear eventIndex startPos eventStart eventEnd isnotDC
                                eventIndex = strcmp(conditionList,uniqueConds(iU)); % create event index
                                startPos = find(eventIndex==1); % create vector of positions in eventTime when event starts
                                eventStart = conditionStartTimes(startPos); % event starts
                                eventEnd = nan(size(startPos)); % preallocate
                                if startPos(end) == length(conditionStartTimes)
                                    eventEnd(1:length(startPos)-1) = conditionStartTimes(startPos(1:length(startPos)-1)+1); % all events but last event end at beginning of next event
                                    eventEnd(length(startPos)) = time(end); % last event ends at end of time for this run
                                else
                                    eventEnd = conditionStartTimes(startPos+1); % each event ends at beginning of next event
                                end
                                isnotDC = logical(nodcTrial(startPos)); % logical index of which trials are not drift corrected

                                % preallocate for speed
                                clear pos 
                                pos.fixTime = cell(length(eventStart));
                                pos.time = cell(length(eventStart));
                                pos.nanTime = cell(length(eventStart));
                                pos.ResFixTime = cell(length(eventStart));
                                pos.saccSTime = cell(length(eventStart));
                                pos.saccETime = cell(length(eventStart));
                                pos.saccBlink = cell(length(eventStart));
                                pos.saccSmall = cell(length(eventStart));
                                pos.saccSaETime = cell(length(eventStart));

                                % find values of interest at each instance of event
                                for iI = 1:length(eventStart); % cycle through instances of event
                                    % raw data stats:
                                    eventDur{iU}(iI) = eventEnd(iI)-eventStart(iI); % calculate duration of each instance of event
                                    pos.time{iI} = find(eventStart(iI)<=time(1:end) & time(1:end)<eventEnd(iI)); % find positions in time during instance of event
                                    pos.nanTime{iI} = find(isnan(dist(pos.time{iI})) == 1); % find positions in where dist == nan during event
                                    nanTime{iU}(iI) = length(pos.nanTime{iI}); % calculate number of timepoints during which there is no data on each instance of event 
                                    pos.fixTime{iI} = find(dist(pos.time{iI}) < fixationRadiusThreshold); % find positions in fixTime during event
                                    fixTime{iU}(iI) = length(pos.fixTime{iI}); % calculate number of timepoints during which subject fixated on each instance of event
                                    eventDistM{iU}(iI) = nanmean(dist(pos.time{iI})); % find mean of distance measures of eye from fixation for timepoints during instance of event
                                    eventDistSD{iU}(iI) = nanstd(dist(pos.time{iI})); % find sd of distance measures of eye from fixation for timepoints during instance of event
                                    eventAngleM{iU}(iI) = nanmean(angle(pos.time{iI})); % find mean of angle measures of eye from fixation for timepoints during instance of event
                                    eventAngleSD{iU}(iI) = nanstd(angle(pos.time{iI})); % find sd of angle measures of eye from fixation for timepoints during instance of event
                                    eventN{iU}(iI) = length(dist(pos.time{iI})); % find number of distance measures for timepoints during instance of event
                                    % drift corrected data stats:
                                    if isnotDC(iI)
                                        pos.ResFixTime{iI} = []; % find positions in time during event where residualDist < fixThreshold
                                        ResFixTime{iU}(iI) = 0; % calculate number of timepoints during which subject fixated on each instance of event
                                        ResDistM{iU}(iI) = nan; % find mean of residuals from linear model for timepoints during instance of event
                                        ResDistSD{iU}(iI) = nan; % find sd of residuals from linear model for timepoints during instance of event
                                        ResAngleM{iU}(iI) = nan; % find mean of residuals from linear model for timepoints during instance of event
                                        ResAngleSD{iU}(iI) = nan; % find sd of residuals from linear model for timepoints during instance of event
                                        ResN{iU}(iI) = 0; % find number of residuals for timepoints during instance of event
                                        nodcTime{iU}(iI) = eventEnd(iI)-eventStart(iI); % amount of time where drift correction was not applied
                                    else % only use stats for drift corrected trials
                                        pos.ResFixTime{iI} = find(abs(dcDist(pos.time{iI})) < fixationRadiusThreshold); % find positions in time during event where residualDist < fixThreshold
                                        ResFixTime{iU}(iI) = length(pos.ResFixTime{iI}); % calculate number of timepoints during which subject fixated on each instance of event
                                        ResDistM{iU}(iI) = nanmean(dcDist(pos.time{iI})); % find mean of residuals from linear model for timepoints during instance of event
                                        ResDistSD{iU}(iI) = nanstd(dcDist(pos.time{iI})); % find sd of residuals from linear model for timepoints during instance of event
                                        ResAngleM{iU}(iI) = nanmean(dcAngle(pos.time{iI})); % find mean of residuals from linear model for timepoints during instance of event
                                        ResAngleSD{iU}(iI) = nanstd(dcAngle(pos.time{iI})); % find sd of residuals from linear model for timepoints during instance of event
                                        ResN{iU}(iI) = length(dcDist(pos.time{iI})); % find number of residuals for timepoints during instance of event
                                        nodcTime{iU}(iI) = 0; % amount of time where drift correction was not applied
                                    end
                                    % saccade and blink stats:
                                    pos.saccSTime{iI} = find(eventStart(iI)<=sacc.Stimes(1:end) & sacc.Stimes(1:end)<eventEnd(iI)); % find positions in sacc.Stimes during event
                                    pos.saccETime{iI} = find(eventStart(iI)<=sacc.Etimes(1:end) & sacc.Etimes(1:end)<eventEnd(iI)); % find positions in sacc.Etimes during event
                                    % add blink count
                                    pos.saccBlink{iI} = pos.saccSTime{iI}(isnan(BSidx(ismember(time,sacc.Stimes(pos.saccSTime{iI}))))); % find positions in sacc.Stimes where there was a blink
                                    blnkN{iU}(iI) = length(pos.saccBlink{iI}); % find number of blinks during instance of event
                                    % remove blinks from saccade count
                                    saccN{iU}(iI) = length(setdiff(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes(pos.saccBlink{iI}))); % find number of saccades (that are not blinks) during instance of event
                                    % create separate saccade count that excludes saccades w/ amp>fixThreshold
                                    pos.saccSmall{iI} = pos.saccSTime{iI}(ismember(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes(sacc.Amps < fixationRadiusThreshold))==1); % find positions in sacc.Stimes where the saccade amplitude < fixThreshold
                                    bigSaccN{iU}(iI) = length(setdiff(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes([pos.saccBlink{iI} pos.saccSmall{iI}])));
                                    pos.saccSaETime{iI} = intersect(pos.saccSTime{iI},pos.saccETime{iI}); % find positions in array of saccades where saccades start and end within instance of event
                                    saccDur{iU}(iI) = sum(sacc.Durs(pos.saccSaETime{iI})); % find the sum of saccade durations which occur within instance of event
                                    if length(pos.saccSTime{iI}) >= length(pos.saccETime{iI}) % find longer pos.sacc*Time array
                                        for st = 1:length(pos.saccSTime{iI}); %cycle through saccades starting during instance of event
                                            % add to saccDur(iI) the within-instance duration of
                                            % saccades starting, but not ending, in event instance 
                                            if st <= length(pos.saccSTime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccSTime{iI}(st)))
                                                saccDur{iU}(iI) = saccDur{iU}(iI) + (eventEnd(iI) - sacc.Stimes(pos.saccSTime{iI}(st)));
                                            % add to saccDur(iI) the within-instance duration of
                                            % saccades ending, but not starting, in event instance
                                            elseif st <= length(pos.saccETime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccETime{iI}(st)))
                                                saccDur{iU}(iI) = saccDur{iU}(iI) + (sacc.Etimes(pos.saccETime{iI}(st)) - eventStart(iI));
                                            end
                                        end
                                    elseif length(pos.saccSTime{iI}) < length(pos.saccETime{iI}) % find longer pos.sacc*Time array
                                        for st = 1:length(pos.saccETime{iI}); %cycle through saccades ending during instance of event
                                            % add to saccDur(i) the within-instance duration of
                                            % saccades starting, but not ending, in event instance 
                                            if st <= length(pos.saccSTime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccSTime{iI}(st)))
                                                saccDur{iU}(iI) = saccDur{iU}(iI) + (eventEnd(iI) - sacc.Stimes(pos.saccSTime{iI}(st)));
                                            % add to saccDur(i) the within-instance duration of
                                            % saccades ending, but not starting, in event instance
                                            elseif st <= length(pos.saccETime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccETime{iI}(st)))
                                                saccDur{iU}(iI) = saccDur{iU}(iI) + (sacc.Etimes(pos.saccETime{iI}(st)) - eventStart(iI));
                                            end
                                        end
                                    end
                                    saccAmpM{iU}(iI) = nanmean(sacc.Amps(pos.saccSTime{iI})); % find mean of amp for timepoints during instance of event
                                    saccAmpSD{iU}(iI) = nanstd(sacc.Amps(pos.saccSTime{iI})); % find sd of amp for timepoints during instance of event
                                    saccAmpN{iU}(iI) = length(sacc.Amps(pos.saccSTime{iI})); % find number of amp measures for timepoints during instance of event
                                    saccVeloM{iU}(iI) = nanmean(sacc.Velo(pos.saccSTime{iI})); % find mean of velocity for timepoints during instance of event
                                    saccVeloSD{iU}(iI) = nanstd(sacc.Velo(pos.saccSTime{iI})); % find sd of velocity for timepoints during instance of event
                                    saccVeloN{iU}(iI) = length(sacc.Velo(pos.saccSTime{iI})); % find number of velocity measures for timepoints during instance of event
                                end

                                % calculate statistics to be saved in data structure
                                if ~isempty(uniqueConds(iU)) % only assign meaningful values to data structure for actual events
                                    % document number of instances of event
                                    eventInstances(iU) = length(eventStart);

                                    % runStats:
                                    % find nan time (%) for event
                                    clear TotalNanTime TotalEventDur
                                    TotalNanTime = nansum(nanTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                    TotalEventDur = nansum(eventDur{iU}(1:end)); % total duration of all instances of event
                                    noTrackTime(iU) = TotalNanTime/TotalEventDur; % document proportion of time fixating across instances of event

                                    % find fixation time (%) for event
                                    clear TotalFixTime
                                    TotalFixTime = nansum(fixTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                    fixationTime(iU) = TotalFixTime/TotalEventDur; % document proportion of time fixating across instances of event

                                    % find distance and angle from fixation stats during event
                                    clear DistGM DistGSD AngleGM AngleGSD
                                    [DistGM, DistGSD] = AK_grandSD(eventDistM{iU},eventDistSD{iU},eventN{iU}); % find mean and SD for all instances of this event
                                    DistM(iU) = DistGM; % document function output
                                    DistSD(iU) = DistGSD; % document function output
                                    [AngleGM, AngleGSD] = AK_grandSD(eventAngleM{iU},eventAngleSD{iU},eventN{iU}); % find mean and SD for all instances of this event
                                    AngleM(iU) = AngleGM; % document function output
                                    AngleSD(iU) = AngleGSD; % document function output

                                    % find drift-corrected fixation time (%) for event
                                    clear TotalResFixTime
                                    TotalResFixTime = nansum(ResFixTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                    driftCorrectedFixationTime(iU) = TotalResFixTime/TotalEventDur; % document proportion of time fixating across instances of event

                                    % find proportion of time for which drift correction was not applied
                                    clear TotalNoDCTime
                                    TotalNoDCTime = nansum(nodcTime{iU}(1:end)); % number of timepoints not drift corrected across all instances of event
                                    nonDriftCorrectedTime(iU) = TotalNoDCTime/TotalEventDur; % document proportion of time fixating across instances of event

                                    % find drift-corrected distance and angle from fixation stats during event
                                    clear ResDistGM ResDistGSD ResAngleGM ResAngleGSD
                                    [ResDistGM, ResDistGSD] = AK_grandSD(ResDistM{iU},ResDistSD{iU},ResN{iU}); % find mean and SD for all instances of this event
                                    driftCorrectedDistM(iU) = ResDistGM; % document function output
                                    driftCorrectedDistSD(iU) = ResDistGSD; % document function output
                                    [ResAngleGM, ResAngleGSD] = AK_grandSD(ResAngleM{iU},ResAngleSD{iU},ResN{iU}); % find mean and SD for all instances of this event
                                    driftCorrectedAngleM(iU) = ResAngleGM; % document function output
                                    driftCorrectedAngleSD(iU) = ResAngleGSD; % document function output

                                    % find saccade stats during event
                                    clear TotalSaccTime
                                    blinkN(iU) = nansum(blnkN{iU}(1:end)); % document sum of number of blinks across instances of event
                                    saccadeN(iU) = nansum(saccN{iU}(1:end)); % document sum of number of sacccades across instances of event
                                    bigSaccadeN(iU) = nansum(bigSaccN{iU}(1:end)); % document sum of number of sacccades outside of fixation threshold across instances of event
                                    TotalSaccTime = nansum(saccDur{iU}(1:end)); % sum saccade durations across instances of event
                                    saccadeTime(iU) = TotalSaccTime/TotalEventDur; % document proportion of time saccading across instances of event
                                    % saccade amplitude stats
                                    clear saccAmpGM saccAmpGSD
                                    [saccAmpGM, saccAmpGSD] = AK_grandSD(saccAmpM{iU},saccAmpSD{iU},saccAmpN{iU}); % find mean and SD for all instances of this event
                                    saccadeAmpM(iU) = saccAmpGM; % document function output
                                    saccadeAmpSD(iU) = saccAmpGSD; % document function output
                                    % saccade velocity stats
                                    clear saccVeloGM saccVeloGSD
                                    [saccVeloGM, saccVeloGSD] = AK_grandSD(saccVeloM{iU},saccVeloSD{iU},saccVeloN{iU}); % find mean and SD for all instances of this event
                                    saccadeVeloM(iU) = saccVeloGM; % document function output
                                    saccadeVeloSD(iU) = saccVeloGSD; % document function output

                                elseif isempty(uniqueConds(iU)) % where uniqueConds(iU) does not exist for file(iA), set stats values = []
                                    eventInstances(iU) = nan;
                                    fixationTime(iU) = nan;
                                    noTrackTime(iU) = nan;
                                    DistM(iU) = nan;
                                    DistSD(iU) = nan;
                                    AngleM(iU) = nan;
                                    AngleSD(iU) = nan;
                                    driftCorrectedFixationTime(iU) = nan;
                                    driftCorrectedDistM(iU) = nan;
                                    driftCorrectedDistSD(iU) = nan;
                                    driftCorrectedAngleM(iU) = nan;
                                    driftCorrectedAngleSD(iU) = nan;
                                    blinkN(iU) = nan;
                                    saccadeN(iU) = nan;
                                    bigSaccadeN(iU) = nan;
                                    saccadeTime(iU) = nan;
                                    saccadeAmpM(iU) = nan;
                                    saccadeAmpSD(iU) = nan;
                                    saccadeVeloM(iU) = nan;
                                    saccadeVeloSD(iU) = nan;
                                end
                            end

                            % calculate statistics across all conditions
                            eventInstances(iU+1) = []; % N/A
                            pos.ALLfixTime = find(dist < fixationRadiusThreshold); % find positions where distance < fixThreshold
                            fixationTime(iU+1) = length(pos.ALLfixTime)/length(time); % calculate fixation time (%) for whole file
                            pos.ALLnanTime = find(isnan(dist) == 1); % find positions in where dist == nan 
                            noTrackTime(iU+1) = length(pos.ALLnanTime)/length(time); % calculate time (%) where there is no tracking
                            DistM(iU+1) = nanmean(dist); % calculate mean distance from fixation
                            DistSD(iU+1) = nanstd(dist); % calculate SD of distance from fixation
                            AngleM(iU+1) = nanmean(angle); % calculate mean angle from fixation
                            AngleSD(iU+1) = nanstd(angle); % calculate SD of angle from fixation
                            pos.ALLResFixTime = find(dcDist < fixationRadiusThreshold); % find positions where drift corrected distance from fixation < fixThreshold
                            driftCorrectedFixationTime(iU+1) = length(pos.ALLResFixTime)/length(time); % calculate drift corrected fixation time (%) for whole file
                            nonDriftCorrectedTime(iU+1) = nansum(nonDriftCorrectedTime(1:iU))/length(time);
                            driftCorrectedDistM(iU+1) = nanmean(dcDist); % calculate mean drift corrected distance from fixation
                            driftCorrectedDistSD(iU+1) = nanstd(dcDist); % calculate SD of drift corrected distance from fixation
                            driftCorrectedAngleM(iU+1) = nanmean(dcAngle); % calculate mean drift corrected distance from fixation
                            driftCorrectedAngleSD(iU+1) = nanstd(dcAngle); % calculate SD of drift corrected distance from fixation
                            blinkN(iU+1) = nansum(blinkN(1:end-1)); % document total number of blinks in whole file
                            saccadeN(iU+1) =  nansum(saccadeN(1:end-1)); % document total number of saccades in whole file
                            bigSaccadeN(iU+1) =  nansum(bigSaccadeN(1:end-1)); % document total number of saccades outside of fixation threshold in whole file
                            saccadeTime(iU+1) = nansum(sacc.Durs)/length(time); % calculate time (%) where subject is in a saccade
                            saccadeAmpM(iU+1) = nanmean(sacc.Amps); % calculate the mean saccade amplitude for the whole file
                            saccadeAmpSD(iU+1) = nanstd(sacc.Amps); % calculate the SD of saccade amplitudes for the whole file
                            saccadeVeloM(iU+1) = nanmean(sacc.Velo); % calculate the mean saccade velocity for the whole file
                            saccadeVeloSD(iU+1) = nanstd(sacc.Velo); % calculate the SD of saccade velocities for the whole file

                            % store stats from asc file iA in data structure
                            subjectData(iSubj).experiment(expIdx).conditionTrialCount{1,thisRun} = eventInstances; % document the number of instances of each condition

                            % designate labels for stats tables
                            clear summaryStatsTable_condLabels statsTable_statLabels
                            runStatsTable_condLabels = [uniqueConds'; {'all'}]';
                            statsTable_statLabels = {'fixationTime','driftCorrectedFixationTime','noDriftCorrectTime','noTrackTime','saccadeTime','distanceMean','distanceSD','angleMean','angleSD','driftCorrectedDistanceMean','driftCorrectedDistanceSD','driftCorrectedAngleMean','driftCorrectedAngleSD','BlinkN','SaccadeN','bigSaccadeN','SaccadeAmplitudeMean','SaccadeAmplitudeSD','SaccadeVelocityMean','SaccadeVelocitySD'};

                            % document table of summary statistics for
                            % each condition in each run
                            subjectData(iSubj).experiment(expIdx).runStats{1,thisRun} = table(fixationTime,driftCorrectedFixationTime,nonDriftCorrectedTime,noTrackTime,saccadeTime,DistM,DistSD,AngleM,AngleSD,driftCorrectedDistM,driftCorrectedDistSD,driftCorrectedAngleM,driftCorrectedAngleSD,blinkN,saccadeN,bigSaccadeN,saccadeAmpM,saccadeAmpSD,saccadeVeloM,saccadeVeloSD,'RowNames',runStatsTable_condLabels,'VariableNames',statsTable_statLabels);

                            % assess data quality for asc file
                            % data is bad quality if all data is missing,
                            % there is no data for at least 40% of the run,
                            % or the standard deviation of the drift
                            % corrected distrance from fixation outside of
                            % saccade events is greater than 1.5 degrees
                            % visual angle
                            if mean(isnan(dist)) == 1 || noTrackTime(end) > .40 || nanstd(dcDist(NoSaccadeTimeIndex)) > 1.5 
                                subjectData(iSubj).experiment(expIdx).quality{1,thisRun} = 0;
                            else
                                subjectData(iSubj).experiment(expIdx).quality{1,thisRun} = 1;
                            end
                        end
                    end
                end
            else
                % fill in fields with nans because there are no asc
                % files for this session directory
                subjectData(iSubj).experiment(expIdx).ascFilename = nan; 
                subjectData(iSubj).experiment(expIdx).behavioralDataFilename = nan; 
                subjectData(iSubj).experiment(expIdx).conditions = nan;
                subjectData(iSubj).experiment(expIdx).conditionTrialCount = nan;
                subjectData(iSubj).experiment(expIdx).runStats = nan;
                subjectData(iSubj).experiment(expIdx).quality = 'could not find asc files';
                disp(['No ' contrast_eye_WildcardStr ' files for subject ' subj_dirs{iSubj}]);
            end
        case 'lz'
            for iSess = 1:length(lz_session_dirs)
                % construct working directory
                use_dir = fullfile(lz_top_dir,subj_dirs{iSubj},lz_session_dirs{iSess});
                % store directory
                subjectData(iSubj).directory{iSess} = use_dir;
                %% deal with motion data:
                % set expIdx and store experiment name
                expIdx = 1;
                subjectData(iSubj).experiment(expIdx).name = 'motion staircase';
                % number of events written to eye tracker during full run
                eventsExpected = 572;
                % set up run counter (runs for a session may span multiple .asc files if the experiment was stopped between runs)
                thisRun = 0;
                % set up mat file number counter (may be distinct from run count if there were bad runs)
                thisMatFileN = 1; % start checking at one
                lastMatFileN = 0; % dummy value
                % generate structure containing information about asc files
                asc_list = dir(fullfile(use_dir,motion_eye_WildcardStr));
                asc_list = AK_sortStruct(asc_list,1,[subj_dirs{iSubj} 'm']); % sort by name
                if ~isempty(asc_list)
                    for iA = 1:length(asc_list); % cycle through asc files
                        % load data
                        clear block mat_name mat_date
                        mat_name = [fullfile(use_dir,asc_list(iA).name(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
                        if ~exist(mat_name,'file') % check whether or not function output has already been saved
                            disp(['parsing ' fullfile(use_dir,asc_list(iA).name)]) % message
                            block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(use_dir,asc_list(iA).name)); % parse asc file
                            save(mat_name,'block'); % save function output
                        else
                            mat_date = dir(mat_name);
                            mat_date = datetime(mat_date.date);
                            if mat_date <= ascParseDate
                                disp(['parsing ' fullfile(use_dir,asc_list(iA).name)]) % message
                                block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(use_dir,asc_list(iA).name)); % parse asc file
                                save(mat_name,'block'); % save function output
                            else
                                disp(['loading ' mat_name]) % message
                                load(mat_name,'block') % load file containing timepts array
                            end
                        end
                        for iB = 1:length(block); % cycle through blocks of trials
                            % only proceed if there are a full set events and if time is continuous
                            if length(block(iB).events(:,1)) == eventsExpected && ~any((cellfun(@str2double,block(iB).timepts(3:end,1))-cellfun(@str2double,block(iB).timepts(2:end-1,1)))~=1)
                                % clear variables
                                clear time Xpixel Ypixel Xdva Ydva dist angle events eventTimes
                                
                                % increase run count
                                thisRun = thisRun + 1;
                                
                                % store asc file name per run per session
                                subjectData(iSubj).experiment(expIdx).ascFilename{iSess,thisRun} = asc_list(iA).name;
                                
                                % extract relevant data from block struct
                                time = cellfun(@str2double,block(iB).timepts(2:end,1)); % set time vector
                                Xpixel = cellfun(@str2double,block(iB).timepts(2:end,2)); % set vector of X position in pixels
                                Ypixel = cellfun(@str2double,block(iB).timepts(2:end,3)); % set vector of Y position in pixels
                                Xdva = (Xpixel-display.center(1))./pixelsPerDegree(1); % set vector of X position as degrees visual angle
                                Ydva = (Ypixel-display.center(2))./pixelsPerDegree(2); % set vector of Y position as degrees visual angle
                                dist = sqrt(Xdva.^2 + Ydva.^2); % calculate vector of distance from fixation as degrees visual angle
                                angle = atan2d(Ydva,Xdva); % calculate vector of polar angles at each timeptpSize = cellfun(@str2double,block(iB).timepts(2:end,4)); % set pupil size vector

                                events = block(iB).events(2:end,2); % store list of events as cell array of strings
                                eventTimes = cellfun(@str2double,block(iB).events(2:end,1)); % set vector of event times
                                events(cellfun(@isempty,events)) = {'no event'}; % fill in empty events
                                
                                % find number of blinks and idx into time array where there are blinks
                                clear blink isBlink
                                blink.N = length(block(iB).blink(:,1))-1;
                                isBlink = zeros(length(time),1);
                                for bl = 2:length(block(iB).blink(:,1)); % cycle through saccades
                                    % get saccade durations, amplitudes, and velocities
                                    blink.Stimes(bl-1) = str2double(block(iB).blink{bl,1});
                                    blink.Etimes(bl-1) = str2double(block(iB).blink{bl,2});
                                    blink.Durs(bl-1) = str2double(block(iB).blink{bl,3});
                                    % find an index in time for each blink
                                    isBlink(blink.Stimes(bl-1)<= time & time<=blink.Etimes(bl-1)) = 1;
                                end

                                % find number of saccades and idx into time array where there are saccades
                                clear sacc isSacc
                                sacc.N = length(block(iB).sacc(:,1))-1;
                                isSacc = zeros(length(time),1);
                                if sacc.N>0 % is there at least one saccade?
                                    for sa = 2:length(block(iB).sacc(:,1)); % cycle through saccades
                                        % get saccade start times, end times, durations, amplitudes, and velocities
                                        sacc.Stimes(sa-1) = str2double(block(iB).sacc{sa,1});
                                        sacc.Etimes(sa-1) = str2double(block(iB).sacc{sa,2});
                                        sacc.Durs(sa-1) = str2double(block(iB).sacc{sa,3}); 
                                        sacc.Amps(sa-1) = str2double(block(iB).sacc{sa,8});
                                        sacc.Velo(sa-1) = str2double(block(iB).sacc{sa,9});
                                        % find an index in time for each saccade
                                        isSacc(sacc.Stimes(sa-1)<= time & time<=sacc.Etimes(sa-1)) = 1;
                                    end
                                else % if there are no saccades
                                    sacc.Stimes = nan;
                                    sacc.Etimes = nan;
                                    sacc.Durs = nan; 
                                    sacc.Amps = nan;
                                    sacc.Velo = nan;
                                end

                                % determine index in time for points which are NOT saccades or blinks
                                clear fullTimeIndex NoSaccadeTimeIndex
                                fullTimeIndex = (1:length(time))';
                                if isfield(sacc,'timeIndex')
                                    NoSaccadeTimeIndex = setdiff(fullTimeIndex,find(isSacc));
                                else
                                    NoSaccadeTimeIndex = fullTimeIndex;
                                end

                                % differentiate saccades with blinks included
                                clear BSidx
                                BSidx = isBlink+isSacc;
                                if length(find(BSidx==0)) > 1
                                    for iSa = 1:sacc.N
                                        if ~isnan(sacc.Stimes(iSa)) && ~isnan(sacc.Etimes(iSa))
                                            clear tempIdxS tempIdxE
                                            tempIdxS = find(time==sacc.Stimes(iSa));
                                            tempIdxE = find(time==sacc.Etimes(iSa));
                                            BSidx(tempIdxS:tempIdxE) = BSidx(tempIdxS:tempIdxE)*max(BSidx(tempIdxS:tempIdxE)); % distinguish blinks from saccades
                                        end
                                    end
                                    BSidx(BSidx>1) = nan; % use to filter blinks from saccade count (eyes open = 0, sacc = 1, blink = nan)
                                end
                                
                                % load behavioral data and match event codes to conditions:
                                clear motionStaircaseN behavioralDataFilename thisRunMatFile conditionIdx conditionList uniqueConds conditionStartTimes condTimeIdx
                                % determine next mat file number, generate
                                % full .mat file name, make sure file exists
                                behavioralDataFilename = [subj_dirs{iSubj} motion_Str num2str(thisMatFileN) '.mat'];
                                thisRunMatFile = fullfile(use_dir,behavioralDataFilename);
                                while exist(thisRunMatFile,'file') == 0 || thisMatFileN == lastMatFileN
                                    thisMatFileN = thisMatFileN + 1;
                                    behavioralDataFilename = [subj_dirs{iSubj} motion_Str num2str(thisMatFileN) '.mat'];
                                    thisRunMatFile = fullfile(use_dir,behavioralDataFilename);
                                end
                                lastMatFileN = thisMatFileN; % in preparation for next run
                                % load
                                load(thisRunMatFile);
                                % match up beginnings of trials (cue times) with condition numbers
                                conditionIdx = runData.conditionOrder; % order of numbered conditions
                                if max(conditionIdx) == 5 % for early version of ssMotion with only 5 conds
                                    conditionIdx(conditionIdx == 5) = 7; % code last condition as catch trials
                                end
                                conditionList = motionCondList(conditionIdx); % list of conditions by name, rather than number
%                                 uniqueConds = unique(conditionList);
                                uniqueConds = motionCondList(ismember(motionCondList,conditionList)); % preserve desired order but use only conditions which were run
                                conditionStartTimes = eventTimes(~cellfun(@isempty,regexp(events,regexptranslate('wildcard','*Trial*_cue_on'))));
                                condTimeIdx = arrayfun(@(x) find(time == x, 1),conditionStartTimes);
                                % store behavioral data filename and run name
                                subjectData(iSubj).experiment(expIdx).behavioralDataFilename{iSess,thisRun} = behavioralDataFilename;
                                subjectData(iSubj).experiment(expIdx).runName{iSess,thisRun} = ['run' num2str(thisRun)];
                                
                                % NEW DRIFT CORRECTION: linear transformations optimized for each trial based on fixation coordinates
                                clear fixations chunkBoundIdx
                                % define what counts as a fixation
                                fixations = AK_defineFixations([Xdva,Ydva,time],fixationRadiusThreshold,fixationMinDur);

                                % separate fixations into chunks by trial and drift correct within each chunk
                                dcX = nan(size(Xdva)); % preallocate variables
                                dcY = nan(size(Ydva));
                                dcIdx = 1;
                                nodcTrial = zeros(length(conditionList),1); % logical to flag trials that receive no drift correction
                                % build chunks as sets of trials that are at least minMillisecondsPerChunk in duration
                                iTrial = 1;
                                chunkIdx = 1;
                                chunkDur = 0;
                                trialChunkIdx = nan(size(conditionStartTimes)); % preallocate as maximum size and trim later
                                while (iTrial + 1) <= length(conditionStartTimes)
                                    % add up trials until duration exceeds minMillisecondsPerChunk
                                    while (iTrial + 1) <= length(conditionStartTimes) && chunkDur < minMillisecondsPerChunk
                                        chunkDur = chunkDur + (conditionStartTimes(iTrial + 1) - conditionStartTimes(iTrial));
                                        iTrial = iTrial + 1;
                                    end
                                    % store trial index of chunk boundaries
                                    trialChunkIdx(chunkIdx) = iTrial;
                                    chunkIdx = chunkIdx + 1; % increment index
                                    chunkDur = 0; % reset chunkDur for next chunk
                                end
                                trialChunkIdx(isnan(trialChunkIdx)) = []; % trim
                                chunkBoundIdx = [condTimeIdx(1);condTimeIdx(trialChunkIdx);length(time)]; % beginning and end indices in time of each chunk
                                for iChunk = 1:length(chunkBoundIdx)-1
                                    %%%test message
                                    disp(['Correcting chunk ' num2str(iChunk) ': timepts ' num2str(chunkBoundIdx(iChunk)) '-' num2str(chunkBoundIdx(iChunk+1))])
                                    %%%
                                    useFix = zeros(length(fixations),1); % reset index of fixations to use for this chunk
                                    chunkFixTime = 0; % reset counter for number of milliseconds spent fixating during this chunk
                                    % determine which fixations to use in the current chunk's drift correction
                                    for iFix = 1:length(fixations)
                                        if isfield(fixations,'timeIdx') && ~isempty(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))) % compare indices for fixation to indices for chunk
                                            useFix(iFix) = 1; % use this fixation for drift correction in the current chunk
                                            chunkFixTime = chunkFixTime+length(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))); % count the amount of time spent fixating within this chunk (assumes 1000hz sample rate)
                                        end
                                    end
                                    useFix = logical(useFix);
                                    % determine whether or not to apply drift correction
                                    clear chunkFix chunkData chunkDCcoords
                                    if any(useFix==1) && chunkFixTime >= 500 % only assign fixation coordinates where there are at least 500ms of fixation in a chunk
                                        % fixation coordinates for this chunk (n fixations by x & y)
                                        chunkFix = vertcat(fixations(useFix).coordinates)'; 
                                        % raw coords for this chunk    
                                        chunkData = horzcat(Xdva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)),Ydva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)))';
                                    else 
                                        % flag non-drift-corrected trials
                                        if iChunk > length(trialChunkIdx) % this last chunk contains all remaining trials in the run
                                            nodcTrial(trialChunkIdx(iChunk-1):length(conditionStartTimes)) = 1;
                                        elseif iChunk == 1 % this first chunk includes all trials from the first until the next chunk boudary
                                            nodcTrial(1:trialChunkIdx(iChunk)) = 1;
                                        else % note all trials in this chunk as normal 
                                            nodcTrial(trialChunkIdx(iChunk-1):trialChunkIdx(iChunk)) = 1;
                                        end
                                        % replace uncorrected data with nans
                                        chunkFix = [];
                                        chunkData = nan(2,length(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)));
                                    end

                                    % apply drift correction to chunk
                                    chunkDCcoords = AK_driftCorrect(chunkData,chunkFix,stimCoords);
                                    % stitch together drift corrected data, removing repeated data points at interior chunkBounds
                                    chunkDCcoords = chunkDCcoords';
                                    if iChunk==1
                                        dcX(dcIdx:length(chunkDCcoords(:,1))) = chunkDCcoords(:,1);
                                        dcY(dcIdx:length(chunkDCcoords(:,2))) = chunkDCcoords(:,2);
                                        dcIdx = dcIdx + length(chunkDCcoords(:,1));
                                    else
                                        dcX(dcIdx:dcIdx + length(chunkDCcoords(2:end,1)) - 1) = chunkDCcoords(2:end,1);
                                        dcY(dcIdx:dcIdx + length(chunkDCcoords(2:end,2)) - 1) = chunkDCcoords(2:end,2);
                                        dcIdx = dcIdx + length(chunkDCcoords(2:end,1));
                                    end
                                end
                                % change to polar coordinates
                                dcDist = sqrt(dcX.^2 + dcY.^2); % calculate vector of distance from fixation as degrees visual angle
                                dcAngle = atan2d(dcY,dcX); % calculate vector of polar angles at each timept
                                
                                if figuresYesNo == 1 % add figures?
                                    % create and name figure
                                    clear figname patchX patchY patchC
                                    figname = strcat('Distance from Fixation Over Time for',' subject:',subj_dirs{iSubj},', session:',lz_session_dirs{iSess},', file:',asc_list(iA).name); % concatenate figure name
                                    figure('Name',figname,'NumberTitle','off'); % generate new figure with proper name

                                    if eventPatchesYesNo == 1 % add patches?
                                        % create patches
                                        patchX = nan(4,length(condTimeIdx)); % preset size of patch input X
                                        patchY = nan(4,length(condTimeIdx)); % preset size of patch input Y
                                        patchC = nan(1,length(condTimeIdx)); % preset size of patch input for color
                                        for iCST = 1:length(conditionStartTimes); % cycle through conditionStartTimes
                                            if iCST ~= length(conditionStartTimes)
                                            patchX(:,iCST) = [conditionStartTimes(iCST) conditionStartTimes(iCST) conditionStartTimes(iCST+1) conditionStartTimes(iCST+1)]'; % set up X coordinates for patches
                                            patchY(:,iCST) = [0 max(dist) max(dist) 0]'; % set up Y coordinates for patches
                                            patchC(1,iCST) = find(strcmp(conditionList(iCST),uniqueConds)==1); % set color index for patches
                                            end
                                        end
                                        
                                        % set color properties of patches
                                        colormap(gray); % set colormap
                                        caxis([1 length(uniqueConds)]); % set number of color intervals

                                        for U = 1:length(uniqueConds) % dumb for loop to make legend work
                                            patch([-100 -100 -99 -99]'+U*iA,[-100 -99 -99 -100]'+U*iA,U)
                                        end

                                        % plot distance from fixation vs time and creat legend
                                        % for stimulus events
                                        patch(patchX,patchY,patchC); % map patches representing events onto axes
                                        legend(uniqueConds);
                                        hold on
                                    end

                                    plot(time,dist,'-b'); % plot distance at time points

                                    % manipulate figure properties
                                    axis_title = strcat(subj_dirs{iSubj},': ',lz_session_dirs{iSess},': ',asc_list(iA).name); % concatenate axis title
                                    title(axis_title)
                                    xlabel('Time (sampled every millisecond)');
                                    ylabel('Distance (degrees visual angle)');
                                    try
                                        axis([min(time(time~=0)) max(time) 0 max(dist)]);
                                    catch
                                        axis([min(time(time~=0)) max(time) 0 5]);
                                    end
                                    set(gca,'XTick',1:50000:max(time))
                                end
                                
                                subjectData(iSubj).experiment(expIdx).conditions{iSess,thisRun} = uniqueConds; % document condition names per asc file

                                % preallocate for speed:
                                % for runStats table
                                clear eventInstances fixationTime noTrackTime DistM DistSD AngleM AngleSD driftCorrectedFixationTime driftCorrectedDistM driftCorrectedDistSD driftCorrectedAngleM driftCorrectedAngleSD saccadeN saccadeTime saccadeAmpM saccadeAmpSD saccadeVeloM saccadeVeloSD
                                eventInstances = zeros(length(uniqueConds)+1,1);
                                fixationTime = zeros(length(uniqueConds)+1,1);
                                noTrackTime = zeros(length(uniqueConds)+1,1);
                                DistM = zeros(length(uniqueConds)+1,1);
                                DistSD = zeros(length(uniqueConds)+1,1);
                                AngleM = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                AngleSD = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                driftCorrectedFixationTime = zeros(length(uniqueConds)+1,1);
                                nonDriftCorrectedTime = zeros(length(uniqueConds)+1,1); % added 9/13/16
                                driftCorrectedDistM = zeros(length(uniqueConds)+1,1);
                                driftCorrectedDistSD = zeros(length(uniqueConds)+1,1);
                                driftCorrectedAngleM = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                driftCorrectedAngleSD = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                blinkN = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                saccadeN = zeros(length(uniqueConds)+1,1);
                                bigSaccadeN = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                saccadeTime = zeros(length(uniqueConds)+1,1);
                                saccadeAmpM = zeros(length(uniqueConds)+1,1);
                                saccadeAmpSD = zeros(length(uniqueConds)+1,1);
                                saccadeVeloM = zeros(length(uniqueConds)+1,1);
                                saccadeVeloSD = zeros(length(uniqueConds)+1,1);
                                % stats by trial
                                clear nanTime fixTime eventDur eventDistM eventDistSD eventAngleM eventAngleSD eventN ResFixTime ResDistM ResDistSD ResAngleM ResAngleSD ResN saccN bigSaccN saccDur saccAmpM saccAmpSD saccAmpN saccVeloM saccVeloSD saccVeloN
                                nanTime = cell(length(uniqueConds),1); 
                                fixTime = cell(length(uniqueConds),1); 
                                eventDur = cell(length(uniqueConds),1); 
                                eventDistM = cell(length(uniqueConds),1); 
                                eventDistSD = cell(length(uniqueConds),1); 
                                eventAngleM = cell(length(uniqueConds),1);  % added 7/27/16
                                eventAngleSD = cell(length(uniqueConds),1);  % added 7/27/16
                                eventN = cell(length(uniqueConds),1); 
                                ResFixTime = cell(length(uniqueConds),1);  
                                ResDistM = cell(length(uniqueConds),1); 
                                ResDistSD = cell(length(uniqueConds),1); 
                                ResAngleM = cell(length(uniqueConds),1);  % added 7/27/16
                                ResAngleSD = cell(length(uniqueConds),1);  % added 7/27/16
                                ResN = cell(length(uniqueConds),1); 
                                nodcTime = cell(length(uniqueConds),1); % added 9/13/16
                                blnkN = cell(length(uniqueConds),1);  % added 7/27/16
                                saccN = cell(length(uniqueConds),1);  
                                bigSaccN = cell(length(uniqueConds),1);  % added 7/27/16
                                saccDur = cell(length(uniqueConds),1);  
                                saccAmpM = cell(length(uniqueConds),1);  
                                saccAmpSD = cell(length(uniqueConds),1); 
                                saccAmpN = cell(length(uniqueConds),1); 
                                saccVeloM = cell(length(uniqueConds),1); 
                                saccVeloSD = cell(length(uniqueConds),1); 
                                saccVeloN = cell(length(uniqueConds),1); 

                                % extract statistics per condition
                                for iU = 1:length(uniqueConds); % cycle through unique events
                                    % find times corresponding to event start and end
                                    clear eventIndex startPos eventStart eventEnd isnotDC
                                    eventIndex = strcmp(conditionList,uniqueConds(iU)); % create event index
                                    startPos = find(eventIndex==1); % create vector of positions in eventTime when event starts
                                    eventStart = conditionStartTimes(startPos); % event starts
                                    eventEnd = nan(size(startPos)); % preallocate
                                    if startPos(end) == length(conditionStartTimes)
                                        eventEnd(1:length(startPos)-1) = conditionStartTimes(startPos(1:length(startPos)-1)+1); % all events but last event end at beginning of next event
                                        eventEnd(length(startPos)) = time(end); % last event ends at end of time for this run
                                    else
                                        eventEnd = conditionStartTimes(startPos+1); % each event ends at beginning of next event
                                    end
                                    isnotDC = logical(nodcTrial(startPos)); % logical index of which trials are not drift corrected

                                    % preallocate for speed
                                    clear pos 
                                    pos.fixTime = cell(length(eventStart));
                                    pos.time = cell(length(eventStart));
                                    pos.nanTime = cell(length(eventStart));
                                    pos.ResFixTime = cell(length(eventStart));
                                    pos.saccSTime = cell(length(eventStart));
                                    pos.saccETime = cell(length(eventStart));
                                    pos.saccBlink = cell(length(eventStart));
                                    pos.saccSmall = cell(length(eventStart));
                                    pos.saccSaETime = cell(length(eventStart));

                                    % find values of interest at each instance of event
                                    for iI = 1:length(eventStart); % cycle through instances of event
                                        % raw data stats:
                                        eventDur{iU}(iI) = eventEnd(iI)-eventStart(iI); % calculate duration of each instance of event
                                        pos.time{iI} = find(eventStart(iI)<=time(1:end) & time(1:end)<eventEnd(iI)); % find positions in time during instance of event
                                        pos.nanTime{iI} = find(isnan(dist(pos.time{iI})) == 1); % find positions in where dist == nan during event
                                        nanTime{iU}(iI) = length(pos.nanTime{iI}); % calculate number of timepoints during which there is no data on each instance of event 
                                        pos.fixTime{iI} = find(dist(pos.time{iI}) < fixationRadiusThreshold); % find positions in fixTime during event
                                        fixTime{iU}(iI) = length(pos.fixTime{iI}); % calculate number of timepoints during which subject fixated on each instance of event
                                        eventDistM{iU}(iI) = nanmean(dist(pos.time{iI})); % find mean of distance measures of eye from fixation for timepoints during instance of event
                                        eventDistSD{iU}(iI) = nanstd(dist(pos.time{iI})); % find sd of distance measures of eye from fixation for timepoints during instance of event
                                        eventAngleM{iU}(iI) = nanmean(angle(pos.time{iI})); % find mean of angle measures of eye from fixation for timepoints during instance of event
                                        eventAngleSD{iU}(iI) = nanstd(angle(pos.time{iI})); % find sd of angle measures of eye from fixation for timepoints during instance of event
                                        eventN{iU}(iI) = length(dist(pos.time{iI})); % find number of distance measures for timepoints during instance of event
                                        % drift corrected data stats:
                                        if isnotDC(iI)
                                            pos.ResFixTime{iI} = []; % find positions in time during event where residualDist < fixThreshold
                                            ResFixTime{iU}(iI) = 0; % calculate number of timepoints during which subject fixated on each instance of event
                                            ResDistM{iU}(iI) = nan; % find mean of residuals from linear model for timepoints during instance of event
                                            ResDistSD{iU}(iI) = nan; % find sd of residuals from linear model for timepoints during instance of event
                                            ResAngleM{iU}(iI) = nan; % find mean of residuals from linear model for timepoints during instance of event
                                            ResAngleSD{iU}(iI) = nan; % find sd of residuals from linear model for timepoints during instance of event
                                            ResN{iU}(iI) = 0; % find number of residuals for timepoints during instance of event
                                            nodcTime{iU}(iI) = eventEnd(iI)-eventStart(iI); % amount of time where drift correction was not applied
                                        else % only use stats for drift corrected trials
                                            pos.ResFixTime{iI} = find(abs(dcDist(pos.time{iI})) < fixationRadiusThreshold); % find positions in time during event where residualDist < fixThreshold
                                            ResFixTime{iU}(iI) = length(pos.ResFixTime{iI}); % calculate number of timepoints during which subject fixated on each instance of event
                                            ResDistM{iU}(iI) = nanmean(dcDist(pos.time{iI})); % find mean of residuals from linear model for timepoints during instance of event
                                            ResDistSD{iU}(iI) = nanstd(dcDist(pos.time{iI})); % find sd of residuals from linear model for timepoints during instance of event
                                            ResAngleM{iU}(iI) = nanmean(dcAngle(pos.time{iI})); % find mean of residuals from linear model for timepoints during instance of event
                                            ResAngleSD{iU}(iI) = nanstd(dcAngle(pos.time{iI})); % find sd of residuals from linear model for timepoints during instance of event
                                            ResN{iU}(iI) = length(dcDist(pos.time{iI})); % find number of residuals for timepoints during instance of event
                                            nodcTime{iU}(iI) = 0; % amount of time where drift correction was not applied
                                        end
                                        % saccade and blink stats:
                                        pos.saccSTime{iI} = find(eventStart(iI)<=sacc.Stimes(1:end) & sacc.Stimes(1:end)<eventEnd(iI)); % find positions in sacc.Stimes during event
                                        pos.saccETime{iI} = find(eventStart(iI)<=sacc.Etimes(1:end) & sacc.Etimes(1:end)<eventEnd(iI)); % find positions in sacc.Etimes during event
                                        % add blink count
                                        pos.saccBlink{iI} = pos.saccSTime{iI}(isnan(BSidx(ismember(time,sacc.Stimes(pos.saccSTime{iI}))))); % find positions in sacc.Stimes where there was a blink
                                        blnkN{iU}(iI) = length(pos.saccBlink{iI}); % find number of blinks during instance of event
                                        % remove blinks from saccade count
                                        saccN{iU}(iI) = length(setdiff(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes(pos.saccBlink{iI}))); % find number of saccades (that are not blinks) during instance of event
                                        % create separate saccade count that excludes saccades w/ amp>fixThreshold
                                        pos.saccSmall{iI} = pos.saccSTime{iI}(ismember(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes(sacc.Amps < fixationRadiusThreshold))==1); % find positions in sacc.Stimes where the saccade amplitude < fixThreshold
                                        bigSaccN{iU}(iI) = length(setdiff(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes([pos.saccBlink{iI} pos.saccSmall{iI}])));
                                        pos.saccSaETime{iI} = intersect(pos.saccSTime{iI},pos.saccETime{iI}); % find positions in array of saccades where saccades start and end within instance of event
                                        saccDur{iU}(iI) = sum(sacc.Durs(pos.saccSaETime{iI})); % find the sum of saccade durations which occur within instance of event
                                        if length(pos.saccSTime{iI}) >= length(pos.saccETime{iI}) % find longer pos.sacc*Time array
                                            for st = 1:length(pos.saccSTime{iI}); %cycle through saccades starting during instance of event
                                                % add to saccDur(iI) the within-instance duration of
                                                % saccades starting, but not ending, in event instance 
                                                if st <= length(pos.saccSTime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccSTime{iI}(st)))
                                                    saccDur{iU}(iI) = saccDur{iU}(iI) + (eventEnd(iI) - sacc.Stimes(pos.saccSTime{iI}(st)));
                                                % add to saccDur(iI) the within-instance duration of
                                                % saccades ending, but not starting, in event instance
                                                elseif st <= length(pos.saccETime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccETime{iI}(st)))
                                                    saccDur{iU}(iI) = saccDur{iU}(iI) + (sacc.Etimes(pos.saccETime{iI}(st)) - eventStart(iI));
                                                end
                                            end
                                        elseif length(pos.saccSTime{iI}) < length(pos.saccETime{iI}) % find longer pos.sacc*Time array
                                            for st = 1:length(pos.saccETime{iI}); %cycle through saccades ending during instance of event
                                                % add to saccDur(i) the within-instance duration of
                                                % saccades starting, but not ending, in event instance 
                                                if st <= length(pos.saccSTime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccSTime{iI}(st)))
                                                    saccDur{iU}(iI) = saccDur{iU}(iI) + (eventEnd(iI) - sacc.Stimes(pos.saccSTime{iI}(st)));
                                                % add to saccDur(i) the within-instance duration of
                                                % saccades ending, but not starting, in event instance
                                                elseif st <= length(pos.saccETime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccETime{iI}(st)))
                                                    saccDur{iU}(iI) = saccDur{iU}(iI) + (sacc.Etimes(pos.saccETime{iI}(st)) - eventStart(iI));
                                                end
                                            end
                                        end
                                        saccAmpM{iU}(iI) = nanmean(sacc.Amps(pos.saccSTime{iI})); % find mean of amp for timepoints during instance of event
                                        saccAmpSD{iU}(iI) = nanstd(sacc.Amps(pos.saccSTime{iI})); % find sd of amp for timepoints during instance of event
                                        saccAmpN{iU}(iI) = length(sacc.Amps(pos.saccSTime{iI})); % find number of amp measures for timepoints during instance of event
                                        saccVeloM{iU}(iI) = nanmean(sacc.Velo(pos.saccSTime{iI})); % find mean of velocity for timepoints during instance of event
                                        saccVeloSD{iU}(iI) = nanstd(sacc.Velo(pos.saccSTime{iI})); % find sd of velocity for timepoints during instance of event
                                        saccVeloN{iU}(iI) = length(sacc.Velo(pos.saccSTime{iI})); % find number of velocity measures for timepoints during instance of event
                                    end

                                    % calculate statistics to be saved in data structure
                                    if ~isempty(uniqueConds(iU)) % only assign meaningful values to data structure for actual events
                                        % document number of instances of event
                                        eventInstances(iU) = length(eventStart);

                                        % runStats:
                                        % find nan time (%) for event
                                        clear TotalNanTime TotalEventDur
                                        TotalNanTime = nansum(nanTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                        TotalEventDur = nansum(eventDur{iU}(1:end)); % total duration of all instances of event
                                        noTrackTime(iU) = TotalNanTime/TotalEventDur; % document proportion of time fixating across instances of event

                                        % find fixation time (%) for event
                                        clear TotalFixTime
                                        TotalFixTime = nansum(fixTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                        fixationTime(iU) = TotalFixTime/TotalEventDur; % document proportion of time fixating across instances of event

                                        % find distance and angle from fixation stats during event
                                        clear DistGM DistGSD AngleGM AngleGSD
                                        [DistGM, DistGSD] = AK_grandSD(eventDistM{iU},eventDistSD{iU},eventN{iU}); % find mean and SD for all instances of this event
                                        DistM(iU) = DistGM; % document function output
                                        DistSD(iU) = DistGSD; % document function output
                                        [AngleGM, AngleGSD] = AK_grandSD(eventAngleM{iU},eventAngleSD{iU},eventN{iU}); % find mean and SD for all instances of this event
                                        AngleM(iU) = AngleGM; % document function output
                                        AngleSD(iU) = AngleGSD; % document function output

                                        % find drift-corrected fixation time (%) for event
                                        clear TotalResFixTime
                                        TotalResFixTime = nansum(ResFixTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                        driftCorrectedFixationTime(iU) = TotalResFixTime/TotalEventDur; % document proportion of time fixating across instances of event

                                        % find proportion of time for which drift correction was not applied
                                        clear TotalNoDCTime
                                        TotalNoDCTime = nansum(nodcTime{iU}(1:end)); % number of timepoints not drift corrected across all instances of event
                                        nonDriftCorrectedTime(iU) = TotalNoDCTime/TotalEventDur; % document proportion of time fixating across instances of event

                                        % find drift-corrected distance and angle from fixation stats during event
                                        clear ResDistGM ResDistGSD ResAngleGM ResAngleGSD
                                        [ResDistGM, ResDistGSD] = AK_grandSD(ResDistM{iU},ResDistSD{iU},ResN{iU}); % find mean and SD for all instances of this event
                                        driftCorrectedDistM(iU) = ResDistGM; % document function output
                                        driftCorrectedDistSD(iU) = ResDistGSD; % document function output
                                        [ResAngleGM, ResAngleGSD] = AK_grandSD(ResAngleM{iU},ResAngleSD{iU},ResN{iU}); % find mean and SD for all instances of this event
                                        driftCorrectedAngleM(iU) = ResAngleGM; % document function output
                                        driftCorrectedAngleSD(iU) = ResAngleGSD; % document function output

                                        % find saccade stats during event
                                        clear TotalSaccTime
                                        blinkN(iU) = nansum(blnkN{iU}(1:end)); % document sum of number of blinks across instances of event
                                        saccadeN(iU) = nansum(saccN{iU}(1:end)); % document sum of number of sacccades across instances of event
                                        bigSaccadeN(iU) = nansum(bigSaccN{iU}(1:end)); % document sum of number of sacccades outside of fixation threshold across instances of event
                                        TotalSaccTime = nansum(saccDur{iU}(1:end)); % sum saccade durations across instances of event
                                        saccadeTime(iU) = TotalSaccTime/TotalEventDur; % document proportion of time saccading across instances of event
                                        % saccade amplitude stats
                                        clear saccAmpGM saccAmpGSD
                                        [saccAmpGM, saccAmpGSD] = AK_grandSD(saccAmpM{iU},saccAmpSD{iU},saccAmpN{iU}); % find mean and SD for all instances of this event
                                        saccadeAmpM(iU) = saccAmpGM; % document function output
                                        saccadeAmpSD(iU) = saccAmpGSD; % document function output
                                        % saccade velocity stats
                                        clear saccVeloGM saccVeloGSD
                                        [saccVeloGM, saccVeloGSD] = AK_grandSD(saccVeloM{iU},saccVeloSD{iU},saccVeloN{iU}); % find mean and SD for all instances of this event
                                        saccadeVeloM(iU) = saccVeloGM; % document function output
                                        saccadeVeloSD(iU) = saccVeloGSD; % document function output

                                    elseif isempty(uniqueConds(iU)) % where uniqueConds(iU) does not exist for file(iA), set stats values = []
                                        eventInstances(iU) = nan;
                                        fixationTime(iU) = nan;
                                        noTrackTime(iU) = nan;
                                        DistM(iU) = nan;
                                        DistSD(iU) = nan;
                                        AngleM(iU) = nan;
                                        AngleSD(iU) = nan;
                                        driftCorrectedFixationTime(iU) = nan;
                                        driftCorrectedDistM(iU) = nan;
                                        driftCorrectedDistSD(iU) = nan;
                                        driftCorrectedAngleM(iU) = nan;
                                        driftCorrectedAngleSD(iU) = nan;
                                        blinkN(iU) = nan;
                                        saccadeN(iU) = nan;
                                        bigSaccadeN(iU) = nan;
                                        saccadeTime(iU) = nan;
                                        saccadeAmpM(iU) = nan;
                                        saccadeAmpSD(iU) = nan;
                                        saccadeVeloM(iU) = nan;
                                        saccadeVeloSD(iU) = nan;
                                    end
                                end

                                % calculate statistics across all conditions
                                eventInstances(iU+1) = []; % N/A
                                pos.ALLfixTime = find(dist < fixationRadiusThreshold); % find positions where distance < fixThreshold
                                fixationTime(iU+1) = length(pos.ALLfixTime)/length(time); % calculate fixation time (%) for whole file
                                pos.ALLnanTime = find(isnan(dist) == 1); % find positions in where dist == nan 
                                noTrackTime(iU+1) = length(pos.ALLnanTime)/length(time); % calculate time (%) where there is no tracking
                                DistM(iU+1) = nanmean(dist); % calculate mean distance from fixation
                                DistSD(iU+1) = nanstd(dist); % calculate SD of distance from fixation
                                AngleM(iU+1) = nanmean(angle); % calculate mean angle from fixation
                                AngleSD(iU+1) = nanstd(angle); % calculate SD of angle from fixation
                                pos.ALLResFixTime = find(dcDist < fixationRadiusThreshold); % find positions where drift corrected distance from fixation < fixThreshold
                                driftCorrectedFixationTime(iU+1) = length(pos.ALLResFixTime)/length(time); % calculate drift corrected fixation time (%) for whole file
                                nonDriftCorrectedTime(iU+1) = nansum(nonDriftCorrectedTime(1:iU))/length(time);
                                driftCorrectedDistM(iU+1) = nanmean(dcDist); % calculate mean drift corrected distance from fixation
                                driftCorrectedDistSD(iU+1) = nanstd(dcDist); % calculate SD of drift corrected distance from fixation
                                driftCorrectedAngleM(iU+1) = nanmean(dcAngle); % calculate mean drift corrected distance from fixation
                                driftCorrectedAngleSD(iU+1) = nanstd(dcAngle); % calculate SD of drift corrected distance from fixation
                                blinkN(iU+1) = nansum(blinkN(1:end-1)); % document total number of blinks in whole file
                                saccadeN(iU+1) =  nansum(saccadeN(1:end-1)); % document total number of saccades in whole file
                                bigSaccadeN(iU+1) =  nansum(bigSaccadeN(1:end-1)); % document total number of saccades outside of fixation threshold in whole file
                                saccadeTime(iU+1) = nansum(sacc.Durs)/length(time); % calculate time (%) where subject is in a saccade
                                saccadeAmpM(iU+1) = nanmean(sacc.Amps); % calculate the mean saccade amplitude for the whole file
                                saccadeAmpSD(iU+1) = nanstd(sacc.Amps); % calculate the SD of saccade amplitudes for the whole file
                                saccadeVeloM(iU+1) = nanmean(sacc.Velo); % calculate the mean saccade velocity for the whole file
                                saccadeVeloSD(iU+1) = nanstd(sacc.Velo); % calculate the SD of saccade velocities for the whole file

                                % store stats from asc file iA in data structure
                                subjectData(iSubj).experiment(expIdx).conditionTrialCount{iSess,thisRun} = eventInstances; % document the number of instances of each condition

                                % designate labels for stats tables
                                clear summaryStatsTable_condLabels statsTable_statLabels
                                runStatsTable_condLabels = [uniqueConds'; {'all'}]';
                                statsTable_statLabels = {'fixationTime','driftCorrectedFixationTime','noDriftCorrectTime','noTrackTime','saccadeTime','distanceMean','distanceSD','angleMean','angleSD','driftCorrectedDistanceMean','driftCorrectedDistanceSD','driftCorrectedAngleMean','driftCorrectedAngleSD','BlinkN','SaccadeN','bigSaccadeN','SaccadeAmplitudeMean','SaccadeAmplitudeSD','SaccadeVelocityMean','SaccadeVelocitySD'};

                                % document table of summary statistics for
                                % each condition in each run
                                subjectData(iSubj).experiment(expIdx).runStats{iSess,thisRun} = table(fixationTime,driftCorrectedFixationTime,nonDriftCorrectedTime,noTrackTime,saccadeTime,DistM,DistSD,AngleM,AngleSD,driftCorrectedDistM,driftCorrectedDistSD,driftCorrectedAngleM,driftCorrectedAngleSD,blinkN,saccadeN,bigSaccadeN,saccadeAmpM,saccadeAmpSD,saccadeVeloM,saccadeVeloSD,'RowNames',runStatsTable_condLabels,'VariableNames',statsTable_statLabels);

                                % assess data quality for asc file
                                % data is bad quality if all data is missing,
                                % there is no data for at least 40% of the run,
                                % or the standard deviation of the drift
                                % corrected distrance from fixation outside of
                                % saccade events is greater than 1.5 degrees
                                % visual angle
                                if mean(isnan(dist)) == 1 || noTrackTime(end) > .40 || nanstd(dcDist(NoSaccadeTimeIndex)) > 1.5 
                                    subjectData(iSubj).experiment(expIdx).quality{iSess,thisRun} = 0;
                                else
                                    subjectData(iSubj).experiment(expIdx).quality{iSess,thisRun} = 1;
                                end                         
                            end
                        end
                    end
                else
                    % fill in fields with nans because there are no asc
                    % files for this session directory
                    subjectData(iSubj).experiment(expIdx).ascFilename{iSess} = nan; 
                    subjectData(iSubj).experiment(expIdx).behavioralDataFilename{iSess} = nan; 
                    subjectData(iSubj).experiment(expIdx).conditions{iSess} = nan;
                    subjectData(iSubj).experiment(expIdx).conditionTrialCount{iSess} = nan;
                    subjectData(iSubj).experiment(expIdx).runStats{iSess} = nan;
                    subjectData(iSubj).experiment(expIdx).quality{iSess} = 'could not find asc files';
                    disp(['No ' motion_eye_WildcardStr ' files for subject ' subj_dirs{iSubj} ', ' lz_session_dirs{iSess}]);
                end
                %% deal with contrast data:
                % set expIdx and store experiment name
                expIdx = 2;
                subjectData(iSubj).experiment(expIdx).name = 'contrast staircase';
                % number of events written to eye tracker during full run
                eventsExpected = 482;
                % set up run counter (runs for a session may span multiple .asc files if the experiment was stopped between runs)
                thisRun = 0;
                % generate structure containing information about asc files
                asc_list = dir(fullfile(use_dir,contrast_eye_WildcardStr)); 
                asc_list = AK_sortStruct(asc_list,1,[subj_dirs{iSubj} 'c']); % sort by name
                if ~isempty(asc_list)
                    for iA = 1:length(asc_list); % cycle through asc files
                        % load data
                        clear block mat_name mat_date
                        mat_name = [fullfile(use_dir,asc_list(iA).name(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
                        if ~exist(mat_name,'file') % check whether or not function output has already been saved
                            disp(['parsing ' fullfile(use_dir,asc_list(iA).name)]) % message
                            block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(use_dir,asc_list(iA).name)); % parse asc file
                            save(mat_name,'block'); % save function output
                        else
                            mat_date = dir(mat_name);
                            mat_date = datetime(mat_date.date);
                            if mat_date <= ascParseDate
                                disp(['parsing ' fullfile(use_dir,asc_list(iA).name)]) % message
                                block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(use_dir,asc_list(iA).name)); % parse asc file
                                save(mat_name,'block'); % save function output
                            else
                                disp(['loading ' mat_name]) % message
                                load(mat_name,'block') % load file containing timepts array
                            end
                        end
                        for iB = 1:length(block); % cycle through blocks of trials
                            % only proceed if there are a full set events and if time is continuous
                            if length(block(iB).events(:,1)) == eventsExpected && ~any((cellfun(@str2double,block(iB).timepts(3:end,1))-cellfun(@str2double,block(iB).timepts(2:end-1,1)))~=1)
                                % clear variables
                                clear time Xpixel Ypixel Xdva Ydva dist angle events eventTimes
                                
                                % increase run count
                                thisRun = thisRun + 1;
                                
                                % store asc file name per run per session
                                subjectData(iSubj).experiment(expIdx).ascFilename{iSess,thisRun} = asc_list(iA).name;
                                
                                % extract relevant data from block struct
                                time = cellfun(@str2double,block(iB).timepts(2:end,1)); % set time vector
                                Xpixel = cellfun(@str2double,block(iB).timepts(2:end,2)); % set vector of X position in pixels
                                Ypixel = cellfun(@str2double,block(iB).timepts(2:end,3)); % set vector of Y position in pixels
                                Xdva = (Xpixel-display.center(1))./pixelsPerDegree(1); % set vector of X position as degrees visual angle
                                Ydva = (Ypixel-display.center(2))./pixelsPerDegree(2); % set vector of Y position as degrees visual angle
                                dist = sqrt(Xdva.^2 + Ydva.^2); % calculate vector of distance from fixation as degrees visual angle
                                angle = atan2d(Ydva,Xdva); % calculate vector of polar angles at each timeptpSize = cellfun(@str2double,block(iB).timepts(2:end,4)); % set pupil size vector

                                events = block(iB).events(2:end,2); % store list of events as cell array of strings
                                eventTimes = cellfun(@str2double,block(iB).events(2:end,1)); % set vector of event times
                                events(cellfun(@isempty,events)) = {'no event'}; % fill in empty events
                                
                                % find number of blinks and idx into time array where there are blinks
                                clear blink isBlink
                                blink.N = length(block(iB).blink(:,1))-1;
                                isBlink = zeros(length(time),1);
                                for bl = 2:length(block(iB).blink(:,1)); % cycle through saccades
                                    % get saccade durations, amplitudes, and velocities
                                    blink.Stimes(bl-1) = str2double(block(iB).blink{bl,1});
                                    blink.Etimes(bl-1) = str2double(block(iB).blink{bl,2});
                                    blink.Durs(bl-1) = str2double(block(iB).blink{bl,3});
                                    % find an index in time for each blink
                                    isBlink(blink.Stimes(bl-1)<= time & time<=blink.Etimes(bl-1)) = 1;
                                end

                                % find number of saccades and idx into time array where there are saccades
                                clear sacc isSacc
                                sacc.N = length(block(iB).sacc(:,1))-1;
                                isSacc = zeros(length(time),1);
                                if sacc.N>0 % is there at least one saccade?
                                    for sa = 2:length(block(iB).sacc(:,1)); % cycle through saccades
                                        % get saccade start times, end times, durations, amplitudes, and velocities
                                        sacc.Stimes(sa-1) = str2double(block(iB).sacc{sa,1});
                                        sacc.Etimes(sa-1) = str2double(block(iB).sacc{sa,2});
                                        sacc.Durs(sa-1) = str2double(block(iB).sacc{sa,3}); 
                                        sacc.Amps(sa-1) = str2double(block(iB).sacc{sa,8});
                                        sacc.Velo(sa-1) = str2double(block(iB).sacc{sa,9});
                                        % find an index in time for each saccade
                                        isSacc(sacc.Stimes(sa-1)<= time & time<=sacc.Etimes(sa-1)) = 1;
                                    end
                                else % if there are no saccades
                                    sacc.Stimes = nan;
                                    sacc.Etimes = nan;
                                    sacc.Durs = nan; 
                                    sacc.Amps = nan;
                                    sacc.Velo = nan;
                                end

                                % determine index in time for points which are NOT saccades or blinks
                                clear fullTimeIndex NoSaccadeTimeIndex
                                fullTimeIndex = (1:length(time))';
                                if isfield(sacc,'timeIndex')
                                    NoSaccadeTimeIndex = setdiff(fullTimeIndex,find(isSacc));
                                else
                                    NoSaccadeTimeIndex = fullTimeIndex;
                                end

                                % differentiate saccades with blinks included
                                clear BSidx
                                BSidx = isBlink+isSacc;
                                if length(find(BSidx==0)) > 1
                                    for iSa = 1:sacc.N
                                        if ~isnan(sacc.Stimes(iSa)) && ~isnan(sacc.Etimes(iSa))
                                            clear tempIdxS tempIdxE
                                            tempIdxS = find(time==sacc.Stimes(iSa));
                                            tempIdxE = find(time==sacc.Etimes(iSa));
                                            BSidx(tempIdxS:tempIdxE) = BSidx(tempIdxS:tempIdxE)*max(BSidx(tempIdxS:tempIdxE)); % distinguish blinks from saccades
                                        end
                                    end
                                    BSidx(BSidx>1) = nan; % use to filter blinks from saccade count (eyes open = 0, sacc = 1, blink = nan)
                                end
                                
                                % load behavioral data and match event codes to conditions:
                                % generate full .mat file name and load
                                clear contrastMatfileN behavioralDataFilename thisRunMatFile conditionIdx conditionList uniqueConds conditionStartTimes condTimeIdx
                                contrastMatfileN = str2double(events{2}(6)); % index into event code (Block#Trial#...) to deal with match eyetracking to behavioral data regardless of counterbalancing
                                behavioralDataFilename = [subj_dirs{iSubj} contrast_lzStr num2str(contrastMatfileN) '_' num2str(iSess) '.mat'];
                                thisRunMatFile = fullfile(use_dir,behavioralDataFilename);
                                load(thisRunMatFile);
                                % match up beginnings of trials (cue times) with condition numbers
                                conditionIdx = runData.conditionOrder; % order of numbered conditions
                                conditionList = contrast_lzCondList(conditionIdx); % list of conditions by name, rather than number
%                                 uniqueConds = unique(conditionList);
                                uniqueConds = contrast_lzCondList(ismember(contrast_lzCondList,conditionList)); % preserve desired order but use only conditions which were run
                                conditionStartTimes = eventTimes(~cellfun(@isempty,regexp(events,regexptranslate('wildcard','*Trial*_cue_on'))));
                                condTimeIdx = arrayfun(@(x) find(time == x, 1),conditionStartTimes);
                                % store behavioral data filename and run name
                                subjectData(iSubj).experiment(expIdx).behavioralDataFilename{iSess,thisRun} = behavioralDataFilename;
                                subjectData(iSubj).experiment(expIdx).runName{iSess,thisRun} = contrast_lzExpList{contrastMatfileN};
                                
                                % NEW DRIFT CORRECTION: linear transformations optimized for each trial based on fixation coordinates
                                clear fixations chunkBoundIdx
                                % define what counts as a fixation
                                fixations = AK_defineFixations([Xdva,Ydva,time],fixationRadiusThreshold,fixationMinDur);

                                % separate fixations into chunks by trial and drift correct within each chunk
                                dcX = nan(size(Xdva)); % preallocate variables
                                dcY = nan(size(Ydva));
                                dcIdx = 1;
                                nodcTrial = zeros(length(conditionList),1); % logical to flag trials that receive no drift correction
                                % build chunks as sets of trials that are at least minMillisecondsPerChunk in duration
                                iTrial = 1;
                                chunkIdx = 1;
                                chunkDur = 0;
                                trialChunkIdx = nan(size(conditionStartTimes)); % preallocate as maximum size and trim later
                                while (iTrial + 1) <= length(conditionStartTimes)
                                    % add up trials until duration exceeds minMillisecondsPerChunk
                                    while (iTrial + 1) <= length(conditionStartTimes) && chunkDur < minMillisecondsPerChunk
                                        chunkDur = chunkDur + (conditionStartTimes(iTrial + 1) - conditionStartTimes(iTrial));
                                        iTrial = iTrial + 1;
                                    end
                                    % store trial index of chunk boundaries
                                    trialChunkIdx(chunkIdx) = iTrial;
                                    chunkIdx = chunkIdx + 1; % increment index
                                    chunkDur = 0; % reset chunkDur for next chunk
                                end
                                trialChunkIdx(isnan(trialChunkIdx)) = []; % trim
                                chunkBoundIdx = [condTimeIdx(1);condTimeIdx(trialChunkIdx);length(time)]; % beginning and end indices in time of each chunk
                                for iChunk = 1:length(chunkBoundIdx)-1
                                    %%%test message
                                    disp(['Correcting chunk ' num2str(iChunk) ': timepts ' num2str(chunkBoundIdx(iChunk)) '-' num2str(chunkBoundIdx(iChunk+1))])
                                    %%%
                                    useFix = zeros(length(fixations),1); % reset index of fixations to use for this chunk
                                    chunkFixTime = 0; % reset counter for number of milliseconds spent fixating during this chunk
                                    % determine which fixations to use in the current chunk's drift correction
                                    for iFix = 1:length(fixations)
                                        if isfield(fixations,'timeIdx') && ~isempty(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))) % compare indices for fixation to indices for chunk
                                            useFix(iFix) = 1; % use this fixation for drift correction in the current chunk
                                            chunkFixTime = chunkFixTime+length(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))); % count the amount of time spent fixating within this chunk (assumes 1000hz sample rate)
                                        end
                                    end
                                    useFix = logical(useFix);
                                    % determine whether or not to apply drift correction
                                    clear chunkFix chunkData chunkDCcoords
                                    if any(useFix==1) && chunkFixTime >= 500 % only assign fixation coordinates where there are at least 500ms of fixation in a chunk
                                        % fixation coordinates for this chunk (n fixations by x & y)
                                        chunkFix = vertcat(fixations(useFix).coordinates)'; 
                                        % raw coords for this chunk    
                                        chunkData = horzcat(Xdva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)),Ydva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)))';
                                    else 
                                        % flag non-drift-corrected trials
                                        if iChunk > length(trialChunkIdx) % this last chunk contains all remaining trials in the run
                                            nodcTrial(trialChunkIdx(iChunk-1):length(conditionStartTimes)) = 1;
                                        elseif iChunk == 1 % this first chunk includes all trials from the first until the next chunk boudary
                                            nodcTrial(1:trialChunkIdx(iChunk)) = 1;
                                        else % note all trials in this chunk as normal 
                                            nodcTrial(trialChunkIdx(iChunk-1):trialChunkIdx(iChunk)) = 1;
                                        end
                                        % replace uncorrected data with nans
                                        chunkFix = [];
                                        chunkData = nan(2,length(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)));
                                    end
                                    % apply drift correction to chunk
                                    chunkDCcoords = AK_driftCorrect(chunkData,chunkFix,stimCoords);
                                    % stitch together drift corrected data, removing repeated data points at interior chunkBounds
                                    chunkDCcoords = chunkDCcoords';
                                    if iChunk==1
                                        dcX(dcIdx:length(chunkDCcoords(:,1))) = chunkDCcoords(:,1);
                                        dcY(dcIdx:length(chunkDCcoords(:,2))) = chunkDCcoords(:,2);
                                        dcIdx = dcIdx + length(chunkDCcoords(:,1));
                                    else
                                        dcX(dcIdx:dcIdx + length(chunkDCcoords(2:end,1)) - 1) = chunkDCcoords(2:end,1);
                                        dcY(dcIdx:dcIdx + length(chunkDCcoords(2:end,2)) - 1) = chunkDCcoords(2:end,2);
                                        dcIdx = dcIdx + length(chunkDCcoords(2:end,1));
                                    end
                                end
                                % change to polar coordinates
                                dcDist = sqrt(dcX.^2 + dcY.^2); % calculate vector of distance from fixation as degrees visual angle
                                dcAngle = atan2d(dcY,dcX); % calculate vector of polar angles at each timept
                                
                                if figuresYesNo == 1 % add figures?
                                    % create and name figure
                                    clear figname patchX patchY patchC
                                    figname = strcat('Distance from Fixation Over Time for',' subject:',subj_dirs{iSubj},', session:',lz_session_dirs{iSess},', file:',asc_list(iA).name); % concatenate figure name
                                    figure('Name',figname,'NumberTitle','off'); % generate new figure with proper name

                                    if eventPatchesYesNo == 1 % add patches?
                                        % create patches
                                        patchX = nan(4,length(condTimeIdx)); % preset size of patch input X
                                        patchY = nan(4,length(condTimeIdx)); % preset size of patch input Y
                                        patchC = nan(1,length(condTimeIdx)); % preset size of patch input for color
                                        for iCST = 1:length(conditionStartTimes); % cycle through conditionStartTimes
                                            if iCST ~= length(conditionStartTimes)
                                            patchX(:,iCST) = [conditionStartTimes(iCST) conditionStartTimes(iCST) conditionStartTimes(iCST+1) conditionStartTimes(iCST+1)]'; % set up X coordinates for patches
                                            patchY(:,iCST) = [0 max(dist) max(dist) 0]'; % set up Y coordinates for patches
                                            patchC(1,iCST) = find(strcmp(conditionList(iCST),uniqueConds)==1); % set color index for patches
                                            end
                                        end
                                        
                                        % set color properties of patches
                                        colormap(gray); % set colormap
                                        caxis([1 length(uniqueConds)]); % set number of color intervals

                                        for U = 1:length(uniqueConds) % dumb for loop to make legend work
                                            patch([-100 -100 -99 -99]'+U*iA,[-100 -99 -99 -100]'+U*iA,U)
                                        end

                                        % plot distance from fixation vs time and creat legend
                                        % for stimulus events
                                        patch(patchX,patchY,patchC); % map patches representing events onto axes
                                        legend(uniqueConds);
                                        hold on
                                    end

                                    plot(time,dist,'-b'); % plot distance at time points

                                    % manipulate figure properties
                                    axis_title = strcat(subj_dirs{iSubj},': ',lz_session_dirs{iSess},': ',asc_list(iA).name); % concatenate axis title
                                    title(axis_title)
                                    xlabel('Time (sampled every millisecond)');
                                    ylabel('Distance (degrees visual angle)');
                                    try
                                        axis([min(time(time~=0)) max(time) 0 max(dist)]);
                                    catch
                                        axis([min(time(time~=0)) max(time) 0 5]);
                                    end
                                    set(gca,'XTick',1:50000:max(time))
                                end
                                
                                subjectData(iSubj).experiment(expIdx).conditions{iSess,thisRun} = uniqueConds; % document condition names per asc file

                                % preallocate for speed:
                                % for runStats table
                                clear eventInstances fixationTime noTrackTime DistM DistSD AngleM AngleSD driftCorrectedFixationTime driftCorrectedDistM driftCorrectedDistSD driftCorrectedAngleM driftCorrectedAngleSD saccadeN saccadeTime saccadeAmpM saccadeAmpSD saccadeVeloM saccadeVeloSD
                                eventInstances = zeros(length(uniqueConds)+1,1);
                                fixationTime = zeros(length(uniqueConds)+1,1);
                                noTrackTime = zeros(length(uniqueConds)+1,1);
                                DistM = zeros(length(uniqueConds)+1,1);
                                DistSD = zeros(length(uniqueConds)+1,1);
                                AngleM = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                AngleSD = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                driftCorrectedFixationTime = zeros(length(uniqueConds)+1,1);
                                nonDriftCorrectedTime = zeros(length(uniqueConds)+1,1); % added 9/13/16
                                driftCorrectedDistM = zeros(length(uniqueConds)+1,1);
                                driftCorrectedDistSD = zeros(length(uniqueConds)+1,1);
                                driftCorrectedAngleM = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                driftCorrectedAngleSD = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                blinkN = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                saccadeN = zeros(length(uniqueConds)+1,1);
                                bigSaccadeN = zeros(length(uniqueConds)+1,1); % added 7/27/16
                                saccadeTime = zeros(length(uniqueConds)+1,1);
                                saccadeAmpM = zeros(length(uniqueConds)+1,1);
                                saccadeAmpSD = zeros(length(uniqueConds)+1,1);
                                saccadeVeloM = zeros(length(uniqueConds)+1,1);
                                saccadeVeloSD = zeros(length(uniqueConds)+1,1);
                                % stats by trial
                                clear nanTime fixTime eventDur eventDistM eventDistSD eventAngleM eventAngleSD eventN ResFixTime ResDistM ResDistSD ResAngleM ResAngleSD ResN saccN bigSaccN saccDur saccAmpM saccAmpSD saccAmpN saccVeloM saccVeloSD saccVeloN
                                nanTime = cell(length(uniqueConds),1); 
                                fixTime = cell(length(uniqueConds),1); 
                                eventDur = cell(length(uniqueConds),1); 
                                eventDistM = cell(length(uniqueConds),1); 
                                eventDistSD = cell(length(uniqueConds),1); 
                                eventAngleM = cell(length(uniqueConds),1);  % added 7/27/16
                                eventAngleSD = cell(length(uniqueConds),1);  % added 7/27/16
                                eventN = cell(length(uniqueConds),1); 
                                ResFixTime = cell(length(uniqueConds),1);  
                                ResDistM = cell(length(uniqueConds),1); 
                                ResDistSD = cell(length(uniqueConds),1); 
                                ResAngleM = cell(length(uniqueConds),1);  % added 7/27/16
                                ResAngleSD = cell(length(uniqueConds),1);  % added 7/27/16
                                ResN = cell(length(uniqueConds),1); 
                                nodcTime = cell(length(uniqueConds),1); % added 9/13/16
                                blnkN = cell(length(uniqueConds),1);  % added 7/27/16
                                saccN = cell(length(uniqueConds),1);  
                                bigSaccN = cell(length(uniqueConds),1);  % added 7/27/16
                                saccDur = cell(length(uniqueConds),1);  
                                saccAmpM = cell(length(uniqueConds),1);  
                                saccAmpSD = cell(length(uniqueConds),1); 
                                saccAmpN = cell(length(uniqueConds),1); 
                                saccVeloM = cell(length(uniqueConds),1); 
                                saccVeloSD = cell(length(uniqueConds),1); 
                                saccVeloN = cell(length(uniqueConds),1); 

                                % extract statistics per condition
                                for iU = 1:length(uniqueConds); % cycle through unique events
                                    % find times corresponding to event start and end
                                    clear eventIndex startPos eventStart eventEnd isnotDC
                                    eventIndex = strcmp(conditionList,uniqueConds(iU)); % create event index
                                    startPos = find(eventIndex==1); % create vector of positions in eventTime when event starts
                                    eventStart = conditionStartTimes(startPos); % event starts
                                    eventEnd = nan(size(startPos)); % preallocate
                                    if startPos(end) == length(conditionStartTimes)
                                        eventEnd(1:length(startPos)-1) = conditionStartTimes(startPos(1:length(startPos)-1)+1); % all events but last event end at beginning of next event
                                        eventEnd(length(startPos)) = time(end); % last event ends at end of time for this run
                                    else
                                        eventEnd = conditionStartTimes(startPos+1); % each event ends at beginning of next event
                                    end
                                    isnotDC = logical(nodcTrial(startPos)); % logical index of which trials are not drift corrected

                                    % preallocate for speed
                                    clear pos 
                                    pos.fixTime = cell(length(eventStart));
                                    pos.time = cell(length(eventStart));
                                    pos.nanTime = cell(length(eventStart));
                                    pos.ResFixTime = cell(length(eventStart));
                                    pos.saccSTime = cell(length(eventStart));
                                    pos.saccETime = cell(length(eventStart));
                                    pos.saccBlink = cell(length(eventStart));
                                    pos.saccSmall = cell(length(eventStart));
                                    pos.saccSaETime = cell(length(eventStart));

                                    % find values of interest at each instance of event
                                    for iI = 1:length(eventStart); % cycle through instances of event
                                        % raw data stats:
                                        eventDur{iU}(iI) = eventEnd(iI)-eventStart(iI); % calculate duration of each instance of event
                                        pos.time{iI} = find(eventStart(iI)<=time(1:end) & time(1:end)<eventEnd(iI)); % find positions in time during instance of event
                                        pos.nanTime{iI} = find(isnan(dist(pos.time{iI})) == 1); % find positions in where dist == nan during event
                                        nanTime{iU}(iI) = length(pos.nanTime{iI}); % calculate number of timepoints during which there is no data on each instance of event 
                                        pos.fixTime{iI} = find(dist(pos.time{iI}) < fixationRadiusThreshold); % find positions in fixTime during event
                                        fixTime{iU}(iI) = length(pos.fixTime{iI}); % calculate number of timepoints during which subject fixated on each instance of event
                                        eventDistM{iU}(iI) = nanmean(dist(pos.time{iI})); % find mean of distance measures of eye from fixation for timepoints during instance of event
                                        eventDistSD{iU}(iI) = nanstd(dist(pos.time{iI})); % find sd of distance measures of eye from fixation for timepoints during instance of event
                                        eventAngleM{iU}(iI) = nanmean(angle(pos.time{iI})); % find mean of angle measures of eye from fixation for timepoints during instance of event
                                        eventAngleSD{iU}(iI) = nanstd(angle(pos.time{iI})); % find sd of angle measures of eye from fixation for timepoints during instance of event
                                        eventN{iU}(iI) = length(dist(pos.time{iI})); % find number of distance measures for timepoints during instance of event
                                        % drift corrected data stats:
                                        if isnotDC(iI)
                                            pos.ResFixTime{iI} = []; % find positions in time during event where residualDist < fixThreshold
                                            ResFixTime{iU}(iI) = 0; % calculate number of timepoints during which subject fixated on each instance of event
                                            ResDistM{iU}(iI) = nan; % find mean of residuals from linear model for timepoints during instance of event
                                            ResDistSD{iU}(iI) = nan; % find sd of residuals from linear model for timepoints during instance of event
                                            ResAngleM{iU}(iI) = nan; % find mean of residuals from linear model for timepoints during instance of event
                                            ResAngleSD{iU}(iI) = nan; % find sd of residuals from linear model for timepoints during instance of event
                                            ResN{iU}(iI) = 0; % find number of residuals for timepoints during instance of event
                                            nodcTime{iU}(iI) = eventEnd(iI)-eventStart(iI); % amount of time where drift correction was not applied
                                        else % only use stats for drift corrected trials
                                            pos.ResFixTime{iI} = find(abs(dcDist(pos.time{iI})) < fixationRadiusThreshold); % find positions in time during event where residualDist < fixThreshold
                                            ResFixTime{iU}(iI) = length(pos.ResFixTime{iI}); % calculate number of timepoints during which subject fixated on each instance of event
                                            ResDistM{iU}(iI) = nanmean(dcDist(pos.time{iI})); % find mean of residuals from linear model for timepoints during instance of event
                                            ResDistSD{iU}(iI) = nanstd(dcDist(pos.time{iI})); % find sd of residuals from linear model for timepoints during instance of event
                                            ResAngleM{iU}(iI) = nanmean(dcAngle(pos.time{iI})); % find mean of residuals from linear model for timepoints during instance of event
                                            ResAngleSD{iU}(iI) = nanstd(dcAngle(pos.time{iI})); % find sd of residuals from linear model for timepoints during instance of event
                                            ResN{iU}(iI) = length(dcDist(pos.time{iI})); % find number of residuals for timepoints during instance of event
                                            nodcTime{iU}(iI) = 0; % amount of time where drift correction was not applied
                                        end
                                        % saccade and blink stats:
                                        pos.saccSTime{iI} = find(eventStart(iI)<=sacc.Stimes(1:end) & sacc.Stimes(1:end)<eventEnd(iI)); % find positions in sacc.Stimes during event
                                        pos.saccETime{iI} = find(eventStart(iI)<=sacc.Etimes(1:end) & sacc.Etimes(1:end)<eventEnd(iI)); % find positions in sacc.Etimes during event
                                        % add blink count
                                        pos.saccBlink{iI} = pos.saccSTime{iI}(isnan(BSidx(ismember(time,sacc.Stimes(pos.saccSTime{iI}))))); % find positions in sacc.Stimes where there was a blink
                                        blnkN{iU}(iI) = length(pos.saccBlink{iI}); % find number of blinks during instance of event
                                        % remove blinks from saccade count
                                        saccN{iU}(iI) = length(setdiff(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes(pos.saccBlink{iI}))); % find number of saccades (that are not blinks) during instance of event
                                        % create separate saccade count that excludes saccades w/ amp>fixThreshold
                                        pos.saccSmall{iI} = pos.saccSTime{iI}(ismember(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes(sacc.Amps < fixationRadiusThreshold))==1); % find positions in sacc.Stimes where the saccade amplitude < fixThreshold
                                        bigSaccN{iU}(iI) = length(setdiff(sacc.Stimes(pos.saccSTime{iI}),sacc.Stimes([pos.saccBlink{iI} pos.saccSmall{iI}])));
                                        pos.saccSaETime{iI} = intersect(pos.saccSTime{iI},pos.saccETime{iI}); % find positions in array of saccades where saccades start and end within instance of event
                                        saccDur{iU}(iI) = sum(sacc.Durs(pos.saccSaETime{iI})); % find the sum of saccade durations which occur within instance of event
                                        if length(pos.saccSTime{iI}) >= length(pos.saccETime{iI}) % find longer pos.sacc*Time array
                                            for st = 1:length(pos.saccSTime{iI}); %cycle through saccades starting during instance of event
                                                % add to saccDur(iI) the within-instance duration of
                                                % saccades starting, but not ending, in event instance 
                                                if st <= length(pos.saccSTime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccSTime{iI}(st)))
                                                    saccDur{iU}(iI) = saccDur{iU}(iI) + (eventEnd(iI) - sacc.Stimes(pos.saccSTime{iI}(st)));
                                                % add to saccDur(iI) the within-instance duration of
                                                % saccades ending, but not starting, in event instance
                                                elseif st <= length(pos.saccETime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccETime{iI}(st)))
                                                    saccDur{iU}(iI) = saccDur{iU}(iI) + (sacc.Etimes(pos.saccETime{iI}(st)) - eventStart(iI));
                                                end
                                            end
                                        elseif length(pos.saccSTime{iI}) < length(pos.saccETime{iI}) % find longer pos.sacc*Time array
                                            for st = 1:length(pos.saccETime{iI}); %cycle through saccades ending during instance of event
                                                % add to saccDur(i) the within-instance duration of
                                                % saccades starting, but not ending, in event instance 
                                                if st <= length(pos.saccSTime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccSTime{iI}(st)))
                                                    saccDur{iU}(iI) = saccDur{iU}(iI) + (eventEnd(iI) - sacc.Stimes(pos.saccSTime{iI}(st)));
                                                % add to saccDur(i) the within-instance duration of
                                                % saccades ending, but not starting, in event instance
                                                elseif st <= length(pos.saccETime{iI}) && isempty(intersect(pos.saccSaETime{iI},pos.saccETime{iI}(st)))
                                                    saccDur{iU}(iI) = saccDur{iU}(iI) + (sacc.Etimes(pos.saccETime{iI}(st)) - eventStart(iI));
                                                end
                                            end
                                        end
                                        saccAmpM{iU}(iI) = nanmean(sacc.Amps(pos.saccSTime{iI})); % find mean of amp for timepoints during instance of event
                                        saccAmpSD{iU}(iI) = nanstd(sacc.Amps(pos.saccSTime{iI})); % find sd of amp for timepoints during instance of event
                                        saccAmpN{iU}(iI) = length(sacc.Amps(pos.saccSTime{iI})); % find number of amp measures for timepoints during instance of event
                                        saccVeloM{iU}(iI) = nanmean(sacc.Velo(pos.saccSTime{iI})); % find mean of velocity for timepoints during instance of event
                                        saccVeloSD{iU}(iI) = nanstd(sacc.Velo(pos.saccSTime{iI})); % find sd of velocity for timepoints during instance of event
                                        saccVeloN{iU}(iI) = length(sacc.Velo(pos.saccSTime{iI})); % find number of velocity measures for timepoints during instance of event
                                    end

                                    % calculate statistics to be saved in data structure
                                    if ~isempty(uniqueConds(iU)) % only assign meaningful values to data structure for actual events
                                        % document number of instances of event
                                        eventInstances(iU) = length(eventStart);

                                        % runStats:
                                        % find nan time (%) for event
                                        clear TotalNanTime TotalEventDur
                                        TotalNanTime = nansum(nanTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                        TotalEventDur = nansum(eventDur{iU}(1:end)); % total duration of all instances of event
                                        noTrackTime(iU) = TotalNanTime/TotalEventDur; % document proportion of time fixating across instances of event

                                        % find fixation time (%) for event
                                        clear TotalFixTime
                                        TotalFixTime = nansum(fixTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                        fixationTime(iU) = TotalFixTime/TotalEventDur; % document proportion of time fixating across instances of event

                                        % find distance and angle from fixation stats during event
                                        clear DistGM DistGSD AngleGM AngleGSD
                                        [DistGM, DistGSD] = AK_grandSD(eventDistM{iU},eventDistSD{iU},eventN{iU}); % find mean and SD for all instances of this event
                                        DistM(iU) = DistGM; % document function output
                                        DistSD(iU) = DistGSD; % document function output
                                        [AngleGM, AngleGSD] = AK_grandSD(eventAngleM{iU},eventAngleSD{iU},eventN{iU}); % find mean and SD for all instances of this event
                                        AngleM(iU) = AngleGM; % document function output
                                        AngleSD(iU) = AngleGSD; % document function output

                                        % find drift-corrected fixation time (%) for event
                                        clear TotalResFixTime
                                        TotalResFixTime = nansum(ResFixTime{iU}(1:end)); % number of timepoints fixated across all instances of event
                                        driftCorrectedFixationTime(iU) = TotalResFixTime/TotalEventDur; % document proportion of time fixating across instances of event

                                        % find proportion of time for which drift correction was not applied
                                        clear TotalNoDCTime
                                        TotalNoDCTime = nansum(nodcTime{iU}(1:end)); % number of timepoints not drift corrected across all instances of event
                                        nonDriftCorrectedTime(iU) = TotalNoDCTime/TotalEventDur; % document proportion of time fixating across instances of event

                                        % find drift-corrected distance and angle from fixation stats during event
                                        clear ResDistGM ResDistGSD ResAngleGM ResAngleGSD
                                        [ResDistGM, ResDistGSD] = AK_grandSD(ResDistM{iU},ResDistSD{iU},ResN{iU}); % find mean and SD for all instances of this event
                                        driftCorrectedDistM(iU) = ResDistGM; % document function output
                                        driftCorrectedDistSD(iU) = ResDistGSD; % document function output
                                        [ResAngleGM, ResAngleGSD] = AK_grandSD(ResAngleM{iU},ResAngleSD{iU},ResN{iU}); % find mean and SD for all instances of this event
                                        driftCorrectedAngleM(iU) = ResAngleGM; % document function output
                                        driftCorrectedAngleSD(iU) = ResAngleGSD; % document function output

                                        % find saccade stats during event
                                        clear TotalSaccTime
                                        blinkN(iU) = nansum(blnkN{iU}(1:end)); % document sum of number of blinks across instances of event
                                        saccadeN(iU) = nansum(saccN{iU}(1:end)); % document sum of number of sacccades across instances of event
                                        bigSaccadeN(iU) = nansum(bigSaccN{iU}(1:end)); % document sum of number of sacccades outside of fixation threshold across instances of event
                                        TotalSaccTime = nansum(saccDur{iU}(1:end)); % sum saccade durations across instances of event
                                        saccadeTime(iU) = TotalSaccTime/TotalEventDur; % document proportion of time saccading across instances of event
                                        % saccade amplitude stats
                                        clear saccAmpGM saccAmpGSD
                                        [saccAmpGM, saccAmpGSD] = AK_grandSD(saccAmpM{iU},saccAmpSD{iU},saccAmpN{iU}); % find mean and SD for all instances of this event
                                        saccadeAmpM(iU) = saccAmpGM; % document function output
                                        saccadeAmpSD(iU) = saccAmpGSD; % document function output
                                        % saccade velocity stats
                                        clear saccVeloGM saccVeloGSD
                                        [saccVeloGM, saccVeloGSD] = AK_grandSD(saccVeloM{iU},saccVeloSD{iU},saccVeloN{iU}); % find mean and SD for all instances of this event
                                        saccadeVeloM(iU) = saccVeloGM; % document function output
                                        saccadeVeloSD(iU) = saccVeloGSD; % document function output

                                    elseif isempty(uniqueConds(iU)) % where uniqueConds(iU) does not exist for file(iA), set stats values = []
                                        eventInstances(iU) = nan;
                                        fixationTime(iU) = nan;
                                        noTrackTime(iU) = nan;
                                        DistM(iU) = nan;
                                        DistSD(iU) = nan;
                                        AngleM(iU) = nan;
                                        AngleSD(iU) = nan;
                                        driftCorrectedFixationTime(iU) = nan;
                                        driftCorrectedDistM(iU) = nan;
                                        driftCorrectedDistSD(iU) = nan;
                                        driftCorrectedAngleM(iU) = nan;
                                        driftCorrectedAngleSD(iU) = nan;
                                        blinkN(iU) = nan;
                                        saccadeN(iU) = nan;
                                        bigSaccadeN(iU) = nan;
                                        saccadeTime(iU) = nan;
                                        saccadeAmpM(iU) = nan;
                                        saccadeAmpSD(iU) = nan;
                                        saccadeVeloM(iU) = nan;
                                        saccadeVeloSD(iU) = nan;
                                    end
                                end

                                % calculate statistics across all conditions
                                eventInstances(iU+1) = []; % N/A
                                pos.ALLfixTime = find(dist < fixationRadiusThreshold); % find positions where distance < fixThreshold
                                fixationTime(iU+1) = length(pos.ALLfixTime)/length(time); % calculate fixation time (%) for whole file
                                pos.ALLnanTime = find(isnan(dist) == 1); % find positions in where dist == nan 
                                noTrackTime(iU+1) = length(pos.ALLnanTime)/length(time); % calculate time (%) where there is no tracking
                                DistM(iU+1) = nanmean(dist); % calculate mean distance from fixation
                                DistSD(iU+1) = nanstd(dist); % calculate SD of distance from fixation
                                AngleM(iU+1) = nanmean(angle); % calculate mean angle from fixation
                                AngleSD(iU+1) = nanstd(angle); % calculate SD of angle from fixation
                                pos.ALLResFixTime = find(dcDist < fixationRadiusThreshold); % find positions where drift corrected distance from fixation < fixThreshold
                                driftCorrectedFixationTime(iU+1) = length(pos.ALLResFixTime)/length(time); % calculate drift corrected fixation time (%) for whole file
                                nonDriftCorrectedTime(iU+1) = nansum(nonDriftCorrectedTime(1:iU))/length(time);
                                driftCorrectedDistM(iU+1) = nanmean(dcDist); % calculate mean drift corrected distance from fixation
                                driftCorrectedDistSD(iU+1) = nanstd(dcDist); % calculate SD of drift corrected distance from fixation
                                driftCorrectedAngleM(iU+1) = nanmean(dcAngle); % calculate mean drift corrected distance from fixation
                                driftCorrectedAngleSD(iU+1) = nanstd(dcAngle); % calculate SD of drift corrected distance from fixation
                                blinkN(iU+1) = nansum(blinkN(1:end-1)); % document total number of blinks in whole file
                                saccadeN(iU+1) =  nansum(saccadeN(1:end-1)); % document total number of saccades in whole file
                                bigSaccadeN(iU+1) =  nansum(bigSaccadeN(1:end-1)); % document total number of saccades outside of fixation threshold in whole file
                                saccadeTime(iU+1) = nansum(sacc.Durs)/length(time); % calculate time (%) where subject is in a saccade
                                saccadeAmpM(iU+1) = nanmean(sacc.Amps); % calculate the mean saccade amplitude for the whole file
                                saccadeAmpSD(iU+1) = nanstd(sacc.Amps); % calculate the SD of saccade amplitudes for the whole file
                                saccadeVeloM(iU+1) = nanmean(sacc.Velo); % calculate the mean saccade velocity for the whole file
                                saccadeVeloSD(iU+1) = nanstd(sacc.Velo); % calculate the SD of saccade velocities for the whole file

                                % store stats from asc file iA in data structure
                                subjectData(iSubj).experiment(expIdx).conditionTrialCount{iSess,thisRun} = eventInstances; % document the number of instances of each condition

                                % designate labels for stats tables
                                clear summaryStatsTable_condLabels statsTable_statLabels
                                runStatsTable_condLabels = [uniqueConds'; {'all'}]';
                                statsTable_statLabels = {'fixationTime','driftCorrectedFixationTime','noDriftCorrectTime','noTrackTime','saccadeTime','distanceMean','distanceSD','angleMean','angleSD','driftCorrectedDistanceMean','driftCorrectedDistanceSD','driftCorrectedAngleMean','driftCorrectedAngleSD','BlinkN','SaccadeN','bigSaccadeN','SaccadeAmplitudeMean','SaccadeAmplitudeSD','SaccadeVelocityMean','SaccadeVelocitySD'};

                                % document table of summary statistics for
                                % each condition in each run
                                subjectData(iSubj).experiment(expIdx).runStats{iSess,thisRun} = table(fixationTime,driftCorrectedFixationTime,nonDriftCorrectedTime,noTrackTime,saccadeTime,DistM,DistSD,AngleM,AngleSD,driftCorrectedDistM,driftCorrectedDistSD,driftCorrectedAngleM,driftCorrectedAngleSD,blinkN,saccadeN,bigSaccadeN,saccadeAmpM,saccadeAmpSD,saccadeVeloM,saccadeVeloSD,'RowNames',runStatsTable_condLabels,'VariableNames',statsTable_statLabels);

                                % assess data quality for asc file
                                % data is bad quality if all data is missing,
                                % there is no data for at least 40% of the run,
                                % or the standard deviation of the drift
                                % corrected distrance from fixation outside of
                                % saccade events is greater than 1.5 degrees
                                % visual angle
                                if mean(isnan(dist)) == 1 || noTrackTime(end) > .40 || nanstd(dcDist(NoSaccadeTimeIndex)) > 1.5 
                                    subjectData(iSubj).experiment(expIdx).quality{iSess,thisRun} = 0;
                                else
                                    subjectData(iSubj).experiment(expIdx).quality{iSess,thisRun} = 1;
                                end
                            end
                        end
                    end
                else
                    % fill in fields with nans because there are no asc
                    % files for this session directory
                    subjectData(iSubj).experiment(expIdx).ascFilename{iSess} = nan; 
                    subjectData(iSubj).experiment(expIdx).behavioralDataFilename{iSess} = nan; 
                    subjectData(iSubj).experiment(expIdx).conditions{iSess} = nan;
                    subjectData(iSubj).experiment(expIdx).conditionTrialCount{iSess} = nan;
                    subjectData(iSubj).experiment(expIdx).runStats{iSess} = nan;
                    subjectData(iSubj).experiment(expIdx).quality{iSess} = 'could not find asc files';
                    disp(['No ' contrast_eye_WildcardStr ' files for subject ' subj_dirs{iSubj} ', ' lz_session_dirs{iSess}]);
                end
            end
    end
end

end

