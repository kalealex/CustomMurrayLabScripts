function [ data ] = AK_GABAinASD_fMRI_EvaluateFixation( subjectCodeCellStr, figuresYesNo, eventPatchesYesNo, saveDataYesNo, ascfileDirectory, homeDirectory)
%AK_GABAinASD_fMRI_EvaluateFixation returns a structure containing summary
%statistics and meta-data per asc file, separating statistics by condition
%if possible.
%   This function will not parse raw edfs. Use the EDF2ASC converter from
%   EyeLink to convert edfs to asc files. 
%   Inputs are as follows:
%    - subjectCodeCellStr: subject code(s) as a cellarray of strings (necessary) 
%    - figuresYesNo: boolean to toggle on or off the generation of
%        figures plotting distance from fixation vs time (optional);
%        defaults to 0
%    - eventPatchesYesNo: boolean to toggle on or off patches representing
%        condition events mapped onto the timecourse plot (optional);
%        defaults to 0
%    - saveDataYesNo: boolean to toggle on or off saving the data
%        structure (optional); defaults to 0
%    - ascfileDirectory: directory of the asc files (optional); defaults to
%        looking in L drive; always assumes subj --> fMRI session -->
%        log_files --> asc files structure
%    - homeDirectory: directory of where to save the data structure
%        (optional); defaults to current directory
%   The output of the function is a data structure containing the fields:
%    - fMRIsession: session names corresponding to rows across each field
%    - ascFilenames: file names corresponding to specific cells across each field
%    - conditions: condition names for each asc file (including 'ModeRecord' and
%        'Pulse' messages as well as 'all' data for the asc file to which each 
%        cell of the field corresponds
%    - conditionInstanceCount: number of instances of each condition
%    - summaryStats: tables of summary statistics for each condition (fixation time,
%        drift-corrected fixation time, no tracking data time, and saccade time as 
%        proportions of the time for each condition; also the mean and SD of distance
%        and polar angle from fixation (drift-corrected and uncorrected), as well as 
%        number of saccades (filtered and unfiltered by amplitude), number of blinks,
%        mean and SD of saccade amplitudes, and mean and SD of saccade velocities).
%    - blockStats: tables of summary statistics for each block of each
%        conditon with the same architecture as the summaryStats table
%    - quality: a description of eye tracking data quality ('none' when
%        data is all nans; 'poor' when less than 60% of the data is usable;
%        'fair' when the SD of drift-corrected distance from fixation is
%        greater than 1 degree visual angle; 'good' otherwise); this is
%        only meant to be descriptive of the extent to which the data is
%        reliable, not the quality of the subject's fixation


% Check arguments
if nargin > 0 && ischar(subjectCodeCellStr)
    subjectCodeCellStr = {subjectCodeCellStr}; % make a single string input into a cell
end
if nargin < 1 || ~iscell(subjectCodeCellStr)
    error('AK_GABAinASD_fMRI_EvaluateFixation requires subject code[s] (cell array of strings) as an argument');
end
if nargin < 6;
   homeDirectory = [];
end
if nargin < 5;
   ascfileDirectory = [];
end
if nargin < 4;
   saveDataYesNo = 0;
end
if nargin < 3;
   eventPatchesYesNo = 0;
end
if nargin < 2;
   figuresYesNo = 0;
end


%% establish directory by subjects and fMRI sessions

if ~exist('ascfileDirectory','var') || isempty(ascfileDirectory);
    ascfileDirectory = 'L:\MurrayLab\ASD\Data';
end
top_dir = ascfileDirectory; % set up base directory
if exist('homeDirectory','var') || ~isempty(homeDirectory);
    home_dir =  homeDirectory; % set directory for saving data
end
subj_dirs = subjectCodeCellStr; % designate folders for invidual participants 
fmri_dirs = {'fMRI1','fMRI2','ftap'}; % designate subfolders
% fmri_dirs = {'ftap'}; % test

%% enter display characteristics for calculating visual angle, pixels per degree, and a fixation threshold

display.dist = 66; % distance to screen in cm
display.width = 33; % width of display in cm
display.height = 24.75; % height of display in cm
display.resolution = [1024 768]; % number of pixels for screen width, height
display.center = display.resolution./2; % screen center in pixel coordinates
display.angles = [2*atand((display.width/2)/display.dist) 2*atand((display.height/2)/display.dist)]; % visual angle of the screen width, height
pixelsPerDegree = display.resolution./display.angles; % pixels per degree visual angle for X dimension and Y dimension

fixationRadiusThreshold = 1; % set threshold for maintaining fixation (degrees visual angle); used to be .5
fixationMinDur = 100; % fixations must be at least this many milliseconds
stimCoords = [0;0]; % coordinates of stimuli on screen (x & y by number of stim)

%% load eyetracking quality table to find data quality information

load(fullfile(top_dir,'EyetrackingQuality.mat'));

%% parse asc files, retrieve timepts cell array, convert data to degrees visual angle, create figures, and calculate descriptive statistics to evaluate fixation

h = waitbar(0,'Creating Directory'); % create waitbar

% preallocate for speed
data(1:length(subj_dirs)) = struct;

for iS = 1:length(subj_dirs); % cycle through subjects
    % skip analysis for disqualified subjects
    if any(strcmp(subj_dirs{iS},dqSubjects))
        disp([dqSubjects{strcmp(subj_dirs{iS},dqSubjects)} 'is disqualified and is not included in analysis']) % message
    else
        for iF = 1:length(fmri_dirs); % cycle through fmri sessions
            clear fmri_dir
            fmri_dir = fullfile(top_dir,subj_dirs{iS},fmri_dirs{iF});

            data(iS).fMRIsession{iF,1} = fmri_dirs{iF}; % document fMRI session name

            if exist(fmri_dir,'dir')==0
                disp(['could not find folder ' fmri_dirs{iF} ' for subject ' subj_dirs{iS}]) % message
            elseif exist(fmri_dir,'dir')==7
                % generate directory in use
                clear use_dir
                use_dir = fullfile(fmri_dir,'log_files'); 
                if exist(use_dir,'dir')==0 % look for directories without log_files subfolder
                    use_dir = fmri_dir; % generate directory in use
                end
                % get list of existing asc files and associated info for this session from quality table
                clear sessQTidx asc_list log_list isData isGoodData
                sessQTidx = strcmp(qualityTable(:,1),subj_dirs{iS}) & strcmp(qualityTable(:,2),fmri_dirs{iF}); % index where quality table rows correspond to subject and session
                asc_list = qualityTable(sessQTidx,4); % list of asc files
                set_list = qualityTable(sessQTidx,3); % list of associated set numbers
                log_list = qualityTable(sessQTidx,5); % list of associated log files
                isData = cell2mat(qualityTable(sessQTidx,6)); % are these asc files from calibrated eyetracking
                isGoodData = cell2mat(qualityTable(sessQTidx,7)); % are these asc files from good quality eyetracking
                if ~all(cellfun(@isempty,asc_list))
                    for iA = 1:length(asc_list); % cycle through asc files
                        % associated file info
                        data(iS).ascFilename{iF,iA} = asc_list{iA}; % document asc file name
                        data(iS).setN{iF,iA} = set_list{iA}; % document associated set number
                        data(iS).logFilename{iF,iA} = log_list{iA}; % document associated logfile name
                        % only analyze real data
                        if isData(iA)
                            
                            wait = iS*iF*iA/(length(subj_dirs)*length(fmri_dirs)*length(asc_list)+1); % set waitbar progress
                            waitbar(wait,h,['processing ' asc_list{iA}(1:end-4)]); % update waitbar

                            clear timepts event mat_name filedate
                            mat_name = [fullfile(use_dir,asc_list{iA}(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
                            % get date modefied for the .mat file
                            filedate = dir(mat_name); 
                            if length(filedate)==1 % there should be one matching file
                                filedate = filedate.date; 
                            else % otherwise reparse
                                filedate = [];
                            end
                            if ~exist(mat_name,'file') % check whether or not function output has already been saved
                                disp(['parsing ' fullfile(use_dir,asc_list{iA})]) % message
                                [timepts, event] = AK_GABAinASD_ascfileParse(fullfile(use_dir,asc_list{iA})); % parse asc file
                                save(mat_name,'timepts','event'); % save function output
                            elseif isempty(filedate) || filedate < datetime(2016,07,29) % parse file if it is older than most recent parsing function
                                disp(['parsing ' fullfile(use_dir,asc_list{iA})]) % message
                                [timepts, event] = AK_GABAinASD_ascfileParse(fullfile(use_dir,asc_list{iA})); % parse asc file
                                save(mat_name,'timepts','event'); % save function output
                            else
                                disp(['loading ' mat_name]) % message
                                load(mat_name,'timepts','event') % load file containing timepts array
                            end

                            % extract and prepare raw data from timepts
                            clear time Xpixel Ypixel Xdva Ydva dist angle events
                            time = cellfun(@str2double,timepts(2:end,1)); % set time vector
                            Xpixel = cellfun(@str2double,timepts(2:end,2)); % set vector of X position in pixels
                            Ypixel = cellfun(@str2double,timepts(2:end,3)); % set vector of Y position in pixels
                            Xdva = (Xpixel-display.center(1))./pixelsPerDegree(1); % set vector of X position as degrees visual angle
                            Ydva = (Ypixel-display.center(2))./pixelsPerDegree(2); % set vector of Y position as degrees visual angle
                            dist = sqrt(Xdva.^2 + Ydva.^2); % calculate vector of distance from fixation as degrees visual angle
                            angle = atan2d(Ydva,Xdva); % calculate vector of polar angles at each timept
                            events = timepts(2:length(timepts),6); % store list of events as cell array of strings

                            % find number of blinks and idx into time array where there are blinks
                            clear blink isBlink
                            blink.N = length(event.blink(:,1))-1;
                            isBlink = zeros(length(time),1);
                            for bl = 2:length(event.blink(:,1)); % cycle through saccades
                                % get saccade durations, amplitudes, and velocities
                                blink.Stimes(bl-1) = str2double(event.blink{bl,1});
                                blink.Etimes(bl-1) = str2double(event.blink{bl,2});
                                blink.Durs(bl-1) = str2double(event.blink{bl,3});
                                % find an index in time for each blink
                                isBlink(blink.Stimes(bl-1)<= time & time<=blink.Etimes(bl-1)) = 1;
                            end

                            % find number of saccades and idx into time array where there are saccades
                            clear sacc isSacc
                            sacc.N = length(event.sacc(:,1))-1;
                            isSacc = zeros(length(time),1);
                            if sacc.N>0 % is there at least one saccade?
                                for sa = 2:length(event.sacc(:,1)); % cycle through saccades
                                    % get saccade start times, end times, durations, amplitudes, and velocities
                                    sacc.Stimes(sa-1) = str2double(event.sacc{sa,1});
                                    sacc.Etimes(sa-1) = str2double(event.sacc{sa,2});
                                    sacc.Durs(sa-1) = str2double(event.sacc{sa,3}); 
                                    sacc.Amps(sa-1) = str2double(event.sacc{sa,8});
                                    sacc.Velo(sa-1) = str2double(event.sacc{sa,9});
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

                            % determine index in time for points which are NOT saccades
                            clear fullTimeIndex NoSaccadeTimeIndex
                            fullTimeIndex = (1:length(time))';
                            if isfield(sacc,'timeIndex')
                                NoSaccadeTimeIndex = setdiff(fullTimeIndex,find(isSacc==1));
                            else
                                NoSaccadeTimeIndex = fullTimeIndex;
                            end

                            % differentiate saccades with blinks included
                            clear BSidx
                            BSidx = isBlink+isSacc;
                            if length(find(BSidx==0)) > 1
                                for iSa = 1:sacc.N
                                    clear tempIdxS tempIdxE
                                    tempIdxS = find(time==sacc.Stimes(iSa)); tempIdxE = find(time==sacc.Etimes(iSa));
                                    BSidx(tempIdxS:tempIdxE) = BSidx(tempIdxS:tempIdxE)*max(BSidx(tempIdxS:tempIdxE));
                                end
                                BSidx(BSidx>1) = nan; % use to filter blinks from saccade count
                            end

                            % create patches and segment information about events
                            clear patchIndex patchEvents uniqueEvents patchTimes patchX patchY patchC
                            patchIndex = find(strcmp(events(1:end-1),events(2:end))==0); % create index of when events change
                            patchEvents = events(patchIndex); % create list of events in order
                            uniqueEvents = unique(patchEvents); % create list of unique events
                            patchTimes = [time(patchIndex); time(end)]; % create vector of time points at which events change
                            patchX = nan(4,length(patchIndex)); % preset size of patch input X
                            patchY = nan(4,length(patchIndex)); % preset size of patch input Y
                            patchC = nan(1,length(patchIndex)); % preset size of patch input for color
                            for iP = 1:length(patchTimes); % cycle through patchTimes
                                if iP ~= length(patchTimes)
                                patchX(:,iP) = [patchTimes(iP) patchTimes(iP) patchTimes(iP+1) patchTimes(iP+1)]'; % set up X coordinates for patches
                                patchY(:,iP) = [0 max(dist) max(dist) 0]'; % set up Y coordinates for patches
                                patchC(1,iP) = find(strcmp(patchEvents(iP),uniqueEvents)==1); % set color index for patches
                                end
                            end


                            % NEW DRIFT CORRECTION: linear transformations optimized for each block based on fixation coordinates
                            clear fixations chunkBoundIdx
                            % define what counts as a fixation
                            fixations = AK_defineFixations([Xdva,Ydva,time],fixationRadiusThreshold,fixationMinDur);

                            % separate fixations into chunks by block and drift correct within each chunk
                            dcX = []; % preallocate variables
                            dcY = [];
                            nodcBlock = zeros(length(patchEvents),1); % logical to flag blocks that receive no drift correction
                            chunkBoundIdx = [1;patchIndex;length(time)]; % beginning and end indices in time of each chunk
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
                                elseif iChunk ~= 1 % otherwise replace uncorrected data with nans
                                    nodcBlock(iChunk-1) = 1; % flag non-drift-corrected blocks
                                    chunkFix = [];
                                    chunkData = nan(2,length(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)));
                                else
                                    chunkFix = [];
                                    chunkData = nan(2,length(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)));
                                end

                                % apply drift correction to chunk
                                chunkDCcoords = AK_driftCorrect(chunkData,chunkFix,stimCoords);
                                % stitch together drift corrected data, removing repeated data points at interior chunkBounds
                                chunkDCcoords = chunkDCcoords';
                                if iChunk==1
                                    dcX = [dcX;chunkDCcoords(:,1)];
                                    dcY = [dcY;chunkDCcoords(:,2)];
                                else
                                    dcX = [dcX;chunkDCcoords(2:end,1)];
                                    dcY = [dcY;chunkDCcoords(2:end,2)];
                                end
                            end
                            % change to polar coordinates
                            dcDist = sqrt(dcX.^2 + dcY.^2); % calculate vector of distance from fixation as degrees visual angle
                            dcAngle = atan2d(dcY,dcX); % calculate vector of polar angles at each timept

                            % OLD DRIFT CORRECTION: linear model to regress out drift; exclude blinks
        %                     clear Xmdl Ymdl residualX residualY dcDist dcAngle
        %                     try
        %                         % run linear model for X and Y dva
        %                         Xmdl = LinearModel.fit(time,Xdva,'linear','Exclude',isnan(Xdva));
        %                         Ymdl = LinearModel.fit(time,Ydva,'linear','Exclude',isnan(Ydva));
        %                         % store vector of residuals
        %                         residualX = Xmdl.Residuals.Raw;
        %                         residualY = Ymdl.Residuals.Raw;
        %                         % convert polar coordinates and use these as drift
        %                         % corrected measures
        %                         dcDist = sqrt(residualX.^2 + residualY.^2); % calculate vector of distance from fixation as degrees visual angle
        %                         dcAngle = atan2d(residualY,residualX); % calculate vector of polar angles at each timept
        %                     catch
        %                         dcDist = nan(length(dist),1);
        %                         dcAngle = nan(length(angle),1);
        %                     end


                            if figuresYesNo == 1 % add figures?
                                % create and name figure
                                clear figname
                                figname = strcat('Distance from Fixation Over Time for',' subject:',subj_dirs{iS},', session:',fmri_dirs{iF},', file:',asc_list{iA}); % concatenate figure name
                                figure('Name',figname,'NumberTitle','off'); % generate new figure with proper name

                                if eventPatchesYesNo == 1 % add patches?
                                    % set color properties of patches
                                    colormap(gray); % set colormap
                                    caxis([1 length(uniqueEvents)]); % set number of color intervals

                                    for U = 1:length(uniqueEvents) % dumb for loop to make legend work
                                        patch([-100 -100 -99 -99]'+U*iA,[-100 -99 -99 -100]'+U*iA,U)
                                    end

                                    % plot distance from fixation vs time and creat legend
                                    % for stimulus events
                                    patch(patchX,patchY,patchC); % map patches representing events onto axes
                                    legend(uniqueEvents);
                                    hold on
                                end

                                plot(time,dist,'-b'); % plot distance at time points

                                % manipulate figure properties
                                axis_title = strcat(subj_dirs{iS},': ',fmri_dirs{iF},': ',asc_list{iA}); % concatenate axis title
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

                            data(iS).conditions{iF,iA} = uniqueEvents; % document condition names per asc file

                            % preallocate for speed:
                            % for summaryStats table
                            clear eventInstances fixationTime noTrackTime DistM DistSD AngleM AngleSD driftCorrectedFixationTime driftCorrectedDistM driftCorrectedDistSD driftCorrectedAngleM driftCorrectedAngleSD saccadeN saccadeTime saccadeAmpM saccadeAmpSD saccadeVeloM saccadeVeloSD
                            eventInstances = zeros(length(uniqueEvents)+1,1);
                            fixationTime = zeros(length(uniqueEvents)+1,1);
                            noTrackTime = zeros(length(uniqueEvents)+1,1);
                            DistM = zeros(length(uniqueEvents)+1,1);
                            DistSD = zeros(length(uniqueEvents)+1,1);
                            AngleM = zeros(length(uniqueEvents)+1,1); % added 7/27/16
                            AngleSD = zeros(length(uniqueEvents)+1,1); % added 7/27/16
                            driftCorrectedFixationTime = zeros(length(uniqueEvents)+1,1);
                            nonDriftCorrectedTime = zeros(length(uniqueEvents)+1,1); % added 9/13/16
                            driftCorrectedDistM = zeros(length(uniqueEvents)+1,1);
                            driftCorrectedDistSD = zeros(length(uniqueEvents)+1,1);
                            driftCorrectedAngleM = zeros(length(uniqueEvents)+1,1); % added 7/27/16
                            driftCorrectedAngleSD = zeros(length(uniqueEvents)+1,1); % added 7/27/16
                            blinkN = zeros(length(uniqueEvents)+1,1); % added 7/27/16
                            saccadeN = zeros(length(uniqueEvents)+1,1);
                            bigSaccadeN = zeros(length(uniqueEvents)+1,1); % added 7/27/16
                            saccadeTime = zeros(length(uniqueEvents)+1,1);
                            saccadeAmpM = zeros(length(uniqueEvents)+1,1);
                            saccadeAmpSD = zeros(length(uniqueEvents)+1,1);
                            saccadeVeloM = zeros(length(uniqueEvents)+1,1);
                            saccadeVeloSD = zeros(length(uniqueEvents)+1,1);
                            % blockStats table
                            clear nanTime fixTime eventDur eventDistM eventDistSD eventAngleM eventAngleSD eventN ResFixTime ResDistM ResDistSD ResAngleM ResAngleSD ResN saccN bigSaccN saccDur saccAmpM saccAmpSD saccAmpN saccVeloM saccVeloSD saccVeloN
                            nanTime = cell(length(uniqueEvents),1); 
                            fixTime = cell(length(uniqueEvents),1); 
                            eventDur = cell(length(uniqueEvents),1); 
                            eventDistM = cell(length(uniqueEvents),1); 
                            eventDistSD = cell(length(uniqueEvents),1); 
                            eventAngleM = cell(length(uniqueEvents),1);  % added 7/27/16
                            eventAngleSD = cell(length(uniqueEvents),1);  % added 7/27/16
                            eventN = cell(length(uniqueEvents),1); 
                            ResFixTime = cell(length(uniqueEvents),1);  
                            ResDistM = cell(length(uniqueEvents),1); 
                            ResDistSD = cell(length(uniqueEvents),1); 
                            ResAngleM = cell(length(uniqueEvents),1);  % added 7/27/16
                            ResAngleSD = cell(length(uniqueEvents),1);  % added 7/27/16
                            ResN = cell(length(uniqueEvents),1); 
                            nodcTime = cell(length(uniqueEvents),1); % added 9/13/16
                            blnkN = cell(length(uniqueEvents),1);  % added 7/27/16
                            saccN = cell(length(uniqueEvents),1);  
                            bigSaccN = cell(length(uniqueEvents),1);  % added 7/27/16
                            saccDur = cell(length(uniqueEvents),1);  
                            saccAmpM = cell(length(uniqueEvents),1);  
                            saccAmpSD = cell(length(uniqueEvents),1); 
                            saccAmpN = cell(length(uniqueEvents),1); 
                            saccVeloM = cell(length(uniqueEvents),1); 
                            saccVeloSD = cell(length(uniqueEvents),1); 
                            saccVeloN = cell(length(uniqueEvents),1); 

                            % extract statistics per condition
                            for iU = 1:length(uniqueEvents); % cycle through unique events
                                % find times corresponding to event start and end
                                clear eventIndex startPos eventStart eventEnd isnotDC
                                eventIndex = strcmp(patchEvents,uniqueEvents(iU)); % create event index
                                startPos = find(eventIndex==1); % create vector of positions in eventTime when event starts
                                eventStart = patchTimes(startPos); % event starts
                                eventEnd = patchTimes(startPos+1); % event ends
                                isnotDC = logical(nodcBlock(startPos)); % logical index of which blocks are not drift corrected

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
                                    else % only use stats for drift corrected blocks
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
                                if ~isempty(uniqueEvents(iU)) % only assign meaningful values to data structure for actual events
                                    % document number of instances of event
                                    eventInstances(iU) = length(eventStart);


                                    % summaryStats:
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


                                    % blockStats:
                                    clear blockFixTime blockResFixTime blockNTT blockSaccTime
                                    blockFixTime = cellfun(@(x,y) x./y,fixTime,eventDur,'UniformOutput',0); % fixation time
                                    blockResFixTime = cellfun(@(x,y) x./y,ResFixTime,eventDur,'UniformOutput',0); % drift-corrected fixation time
                                    blockNTT = cellfun(@(x,y) x./y,nanTime,eventDur,'UniformOutput',0); % no track time
                                    blockNDCT = cellfun(@(x,y) x./y,nodcTime,eventDur,'UniformOutput',0); % non drift corrected time
                                    blockSaccTime = cellfun(@(x,y) x./y,saccDur,eventDur,'UniformOutput',0); % saccade time

                                elseif isempty(uniqueEvents(iU)) % where uniqueEvents(iU) does not exist for file(iA), set stats values = []
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

                                    blockFixTime(iU) = {nan};
                                    blockResFixTime(iU) = {nan};
                                    blockNTT(iU) = {nan};
                                    blockSaccTime(iU) = {nan};
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
                            data(iS).conditionBlockCount{iF,iA} = eventInstances; % document the number of instances of each condition

                            % designate labels for stats tables
                            clear summaryStatsTable_condLabels statsTable_statLabels
                            summaryStatsTable_condLabels = [uniqueEvents; {'all'}]';
                            blockStatsTable_condLabels = uniqueEvents';
                            statsTable_statLabels = {'fixationTime','driftCorrectedFixationTime','noDriftCorrectTime','noTrackTime','saccadeTime','distanceMean','distanceSD','angleMean','angleSD','driftCorrectedDistanceMean','driftCorrectedDistanceSD','driftCorrectedAngleMean','driftCorrectedAngleSD','BlinkN','SaccadeN','bigSaccadeN','SaccadeAmplitudeMean','SaccadeAmplitudeSD','SaccadeVelocityMean','SaccadeVelocitySD'};

                            % document table of summary statistics for each condition
                            data(iS).summaryStats{iF,iA} = table(fixationTime,driftCorrectedFixationTime,nonDriftCorrectedTime,noTrackTime,saccadeTime,DistM,DistSD,AngleM,AngleSD,driftCorrectedDistM,driftCorrectedDistSD,driftCorrectedAngleM,driftCorrectedAngleSD,blinkN,saccadeN,bigSaccadeN,saccadeAmpM,saccadeAmpSD,saccadeVeloM,saccadeVeloSD,'RowNames',summaryStatsTable_condLabels,'VariableNames',statsTable_statLabels);

                            % document table of summary statistics for each block by condition
                            data(iS).blockStats{iF,iA} = table(blockFixTime,blockResFixTime,blockNDCT,blockNTT,blockSaccTime,eventDistM,eventDistSD,eventAngleM,eventAngleSD,ResDistM,ResDistSD,ResAngleM,ResAngleSD,blnkN,saccN,bigSaccN,saccAmpM,saccAmpSD,saccVeloM,saccVeloSD,'RowNames',blockStatsTable_condLabels,'VariableNames',statsTable_statLabels);                    

                            % assess data quality for asc file
                            data(iS).quality{iF,iA} = isGoodData(iA);
%                             if mean(isnan(dist)) == 1; % is there actual data?
%                                 data(iS).quality{iF,iA} = 'none';
%                             elseif noTrackTime(end) > .40; % is less than 60% of the file actual data?
%                                 data(iS).quality{iF,iA} = 'poor';
%                             elseif nanstd(dcDist(NoSaccadeTimeIndex)) > 1; % is the SD of drift-corrected distance from fixation, excluding saccades, is greater than 1 degree visual angle
%                                 data(iS).quality{iF,iA} = 'fair';
%                             else 
%                                 data(iS).quality{iF,iA} = 'good';
%                             end
                        else
                            data(iS).conditions{iF,iA} = nan;
                            data(iS).conditionBlockCount{iF,iA} = nan;
                            data(iS).summaryStats{iF,iA} = nan;
                            data(iS).blockStats{iF,iA} = nan;
                            data(iS).quality{iF,iA} = 0;
                        end
                    end
                else
                    data(iS).quality{iF} = 'could not find asc files';
                    data(iS).ascFilename{iF} = nan; % document asc file name
                    data(iS).setN{iF} = nan; % document associated set number
                    data(iS).logFilename{iF} = nan; % document associated logfile name
                    disp(['could not find .asc files for ' fmri_dirs{iF} ' for subject ' subj_dirs{iS}]) % message
                end
            end
        end
    end
end

% close waitbar
wait = 1; % set waitbar progress
waitbar(wait,h,'Finishing'); % update waitbar
% save data structure
if saveDataYesNo == 1 && (exist('home_dir','var') || ~isempty(home_dir));
    data_name = fullfile(home_dir,'EyeTrackingDataStructure.mat'); % create string to name .mat file which saves data
    disp(['saving ' data_name]) % message
    save(data_name,'data'); % save stats
end
pause(.25)
close(h);


end

