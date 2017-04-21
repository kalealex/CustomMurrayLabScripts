function AK_makeEyetrackingMovie( edfName, directory, isPsychophysics, settings, saveParams )
%AK_makeEyetrackingMovie makes a movie out of eye tracking data by cutting
%the time series data into chunks of at least 10 seconds, based on trials
%in psychophysics or blocks in fMRI, for the prupose of drift correcting
%and plotted both raw and drift corrected data. The movie consists of many
%parts saved in a specified directory. The movie consists of many frames in
%which gaze coordinates are plotted on a Cartesian plane with tick marks
%indicating degrees visual angle. The gaze coordinates are marked bya pair
%of color coded dots, where the white dot is raw (non drift corrected)
%data, and the colored dot is drift corrected (green = fixation; red =
%saccade; blue = blink; black = default). Condition, time (seconds), and
%counts of blinks, small saccades, big saccades, and breaks in fixation are
%all noted in text boxes on the screen. See details below for specific
%details of use.
%   INPUT:
%       edfName: the name of the eyetracking data file as a string (without
%           the file extension); note that the function processes asc
%           files, not edfs, so BEFORE RUNNING, CHECK THAT THE EDF HAS BEEN
%           CONVERTED TO ASC FORMAT USING EDF2ASC (an application from SR
%           Research)
%       directory: the directoy the file is in as a string 
%       isPsychophysics: a boolean (true if the eyetracking comes from
%           psychophysics at CHDD; false if eyetracking comes from fMRI)
%       settings[optional]: a structure with the following fields:
%           driftCorrect: a boolean (do drift correction?); defaults to
%               true
%           fixationRadiusThreshold: a criterion used to define the radius
%               of a cluster of gaze coordinates to be consisdered a
%               fixation; this is used in identifying fixations for drift
%               correction, identifying breaks in fixation, and separating
%               saccades into small (< fixationRadiusThreshold) and big (>
%               fixationRadiusThreshold); defaults to 1 degree visual angle
%           fixationMinDur: a criterion used to define the minimum duration
%               of fixations in milliseconds; this is used only for drift
%               correction; defaults to 100 ms
%           stimCoords: coordinates of stimuli on screen (x & y by number
%               of stim); defaults to [0 0] (we expect them to fixate)
%           axesLimit: boundaries of the movie display will be plus and
%               minus this value, in degrees visual angle, from the center
%               of the screen; defaults to 5 degrees visual angle
%       saveParams[optional: a structure with the following fields:
%           directory: the directory to save the movie parts in; defaults
%               to [directory(second argument) '/' edfName(first argument)
%               '_EyetrackingMovie'] and will create this folder if it
%               doesn't exist
%           frameRate: the frame rate of the movie to be created; defaults
%               to 60 hz
%           quality: the desired percentage of possible resolution in the
%               saved movie; defaults to 50%, which is sufficient 

%% check input

% file and location selection
if nargin < 3
    error(['AK_makeEyetrackingMovie requires as input at least three arguments: (1) the name of the edf you wish to make a movie of as a string;'...
        '(2) the directory of the edf as a string; and (3) a boolean (true if edf was produced during psychophysics on the CHDD setup,'...
        'false if the edf was produced during fMRI)'])
end
% settings
if nargin < 4
    % drift correction params
    settings.driftCorrect = true; % logical (do drift correction?)
    settings.fixationRadiusThreshold = 1; % used to be .5
    settings.fixationMinDur = 100; % fixations must be at least this many milliseconds
    settings.stimCoords = [0;0]; % coordinates of stimuli on screen (x & y by number of stim)
    % axis limits
    settings.axesLimit = 5; % display boundaries as distance from fixation (dva)
end
% movie saving params
if nargin < 5
    % designate save dir
    folderName = [edfName '_EyetrackingMovie'];
    saveParams.directory = fullfile(directory,folderName);
    % set other save params
    saveParams.frameRate = 60; % save function defaults to 30 hz otherwise
    saveParams.quality = 50; % save function defaults to 75 otherwise
end
% make sure directory to save movie exists
if exist(saveParams.directory,'file') == 0
    mkdir(saveParams.directory);
end

%% enter display characteristics for calculating visual angle, pixels per degree

if ~isPsychophysics
    % for fMRI
    display.dist = 66; % distance to screen in cm
    display.width = 33; % width of display in cm
    display.height = 24.75; % height of display in cm
    display.resolution = [1024 768]; % number of pixels for screen width, height
    display.center = display.resolution./2; % screen center in pixel coordinates
    display.angles = [2*atand((display.width/2)/display.dist) 2*atand((display.height/2)/display.dist)]; % visual angle of the screen width, height
    pixelsPerDegree = display.resolution./display.angles; % pixels per degree visual angle for X dimension and Y dimension
else
    % for CHDD Psychophysics setup
    display.dist = 66; % distance to screen in cm
    display.width = 36.5; % width of display in cm
    display.height = 27; % height of display in cm
    display.resolution = [799 599]; % number of pixels for screen width, height
    display.center = display.resolution./2; % screen center in pixel coordinates
    display.angles = [2*atand((display.width/2)/display.dist) 2*atand((display.height/2)/display.dist)]; % visual angle of the screen width, height
    pixelsPerDegree = display.resolution./display.angles; % pixels per degree visual angle for X dimension and Y dimension
end

%% make directory and load file

% generate file names
mat_filename = fullfile(directory,[edfName '.mat']);
asc_filename = fullfile(directory,[edfName '.asc']);

% everything from this point is slightly different for fMRI than it is for
% psychophysics
if isPsychophysics
    % load raw data
    if ~exist(mat_filename,'file') % check whether or not function output has already been saved
        disp(['parsing ' asc_filename]) % message
        block = AK_GABAinASD_ascfileParse(asc_filename); % parse asc file
        save(mat_filename,'block'); % save function output
    else
        disp(['loading ' mat_filename]) % message
        load(mat_filename,'block') % load file containing timepts array
    end
    %% processing: create indices flagging fixation, saccades, and blinks; create indices marking conditions and extract list of condition names; drift correct
    
    % preallocate
    runIdx = 0;
    % cycle through blocks ignoring practice trials
    for iB = 1:length(block)
        if ~any(~cellfun(@isempty,regexp(block(iB).events(:,2),regexptranslate('wildcard','*Prac*'),'once'))) % if there are no events matching the wildcard '*Prac*'
            % increment run counter
            runIdx = runIdx + 1;
            
            % extract and prepare raw data from timepts
            time = cellfun(@str2double,block(iB).timepts(2:end,1)); % set time vector
            Xpixel = cellfun(@str2double,block(iB).timepts(2:end,2)); % set vector of X position in pixels
            Ypixel = cellfun(@str2double,block(iB).timepts(2:end,3)); % set vector of Y position in pixels
            Xdva = (Xpixel-display.center(1))./pixelsPerDegree(1); % set vector of X position as degrees visual angle
            Ydva = (Ypixel-display.center(2))./pixelsPerDegree(2); % set vector of Y position as degrees visual angle
            events = block(iB).events(2:length(block(iB).events),2); % store list of events as cell array of strings
            eventTimes = cellfun(@str2double,block(iB).events(2:length(block(iB).events),1)); % store list of times associated with events as a vector
            
            % blinks and saccades:
            % find number of blinks and idx into time array where there are blinks
            blink.N = length(block(iB).blink(:,1))-1;
            isBlink = zeros(length(time),1);
            if blink.N>0 % is there at least one blink?
                for bl = 2:length(block(iB).blink(:,1)); % cycle through saccades
                    % get blink start times and end times
                    blink.Stimes(bl-1) = str2double(block(iB).blink{bl,1});
                    blink.Etimes(bl-1) = str2double(block(iB).blink{bl,2});
                    % find an index in time for each blink
                    isBlink(blink.Stimes(bl-1)<= time & time<=blink.Etimes(bl-1)) = 1;
                end
            else % if there are no blinks
                blink.Stimes = nan;
                blink.Etimes = nan;
            end
            % find number of saccades and idx into time array where there are saccades
            sacc.N = length(block(iB).sacc(:,1))-1;
            isSacc = zeros(length(time),1);
            if sacc.N>0 % is there at least one saccade?
                for sa = 2:length(block(iB).sacc(:,1)); % cycle through saccades
                    % get saccade start times and end times
                    sacc.Stimes(sa-1) = str2double(block(iB).sacc{sa,1});
                    sacc.Etimes(sa-1) = str2double(block(iB).sacc{sa,2});
                    sacc.Amps(sa-1) = str2double(block(iB).sacc{sa,8});
                    % find an index in time for each saccade
                    if sacc.Amps(sa-1) > settings.fixationRadiusThreshold % for big saccades
                        isSacc(sacc.Stimes(sa-1)<= time & time<=sacc.Etimes(sa-1)) = 1;
                    else % for small saccades
                        isSacc(sacc.Stimes(sa-1)<= time & time<=sacc.Etimes(sa-1)) = .75;
                    end
                end
            else % if there are no saccades
                sacc.Stimes = nan;
                sacc.Etimes = nan;
            end
            % differentiate small saccades, big saccades and blinks
            BSidx = isBlink+isSacc;
            if length(find(BSidx==0)) > 1
                for iSa = 1:sacc.N
                    clear tempIdxS tempIdxE
                    tempIdxS = find(time==sacc.Stimes(iSa)); tempIdxE = find(time==sacc.Etimes(iSa));
                    BSidx(tempIdxS:tempIdxE) = BSidx(tempIdxS:tempIdxE)*max(BSidx(tempIdxS:tempIdxE));
                end
                BSidx(BSidx>1) = nan; % use to distinguish blinks from saccades
            end

            % segment information about events and create indices
            clear eventTrialIdx allEvents trialStartTimes trialTimeIdx
            eventTrialIdx = ~cellfun(@isempty,regexp(events,regexptranslate('wildcard','*Trial*_cue_on')));
            allEvents = events(eventTrialIdx);
            trialStartTimes = eventTimes(eventTrialIdx);
            trialTimeIdx = arrayfun(@(x) find(time == x, 1),trialStartTimes);

            % drift correction of coordinates and fixation index
            % preallocate
            fixIdx = zeros(length(time),1);
            if settings.driftCorrect
                % NEW DRIFT CORRECTION: linear transformations optimized for each block based on fixation coordinates
                % define what counts as a fixation
                fixations = AK_defineFixations([Xdva,Ydva,time],settings.fixationRadiusThreshold,settings.fixationMinDur);

                % separate fixations into chunks by block and drift correct within each chunk
                dcX = []; % preallocate variables
                dcY = [];
                dcIdx = 1;
                % build chunks as sets of trials that are at least minMillisecondsPerChunk in duration
                iTrial = 1;
                chunkIdx = 1;
                chunkDur = 0;
                trialChunkIdx = nan(size(trialStartTimes)); % preallocate as maximum size and trim later
                while (iTrial + 1) <= length(trialStartTimes)
                    % add up trials until duration exceeds minMillisecondsPerChunk
                    while (iTrial + 1) <= length(trialStartTimes) && chunkDur < 10006
                        chunkDur = chunkDur + (trialStartTimes(iTrial + 1) - trialStartTimes(iTrial));
                        iTrial = iTrial + 1;
                    end
                    % store trial index of chunk boundaries
                    trialChunkIdx(chunkIdx) = iTrial;
                    chunkIdx = chunkIdx + 1; % increment index
                    chunkDur = 0; % reset chunkDur for next chunk
                end
                trialChunkIdx(isnan(trialChunkIdx)) = []; % trim
                chunkBoundIdx = [trialTimeIdx(1);trialTimeIdx(trialChunkIdx);length(time)]; % beginning and end indices in time of each chunk
                for iChunk = 1:length(chunkBoundIdx)-1
                    useFix = zeros(length(fixations),1); % reset index of fixations to use for this chunk
                    % determine which fixations to use in the current chunk's drift correction
                    for iFix = 1:length(fixations)
                        if ~isempty(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))) % compare indices for fixation to indices for chunk
                            useFix(iFix) = 1; % use this fixation for drift correction in the current chunk
                        end    
                    end
                    useFix = logical(useFix);
                    % prepare drift correction function input
                    clear chunkFix chunkData chunkDCcoords
                    chunkFix = vertcat(fixations(useFix).coordinates)'; % fixation coordinates for this chunk (n fixations by x & y)
                    chunkData = horzcat(Xdva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)),Ydva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)))';

                    % apply drift correction to chunk
                    chunkDCcoords = AK_driftCorrect(chunkData,chunkFix,settings.stimCoords);
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

                % create fixation idx (zero is fixation)
                fixIdx(sqrt(dcX.^2 + dcY.^2) > settings.fixationRadiusThreshold) = 1;
            else
                % create fixation idx (zero is fixation)
                fixIdx(sqrt(Xdva.^2 + Ydva.^2) > settings.fixationRadiusThreshold) = 1;
            end

            % trim data prior to actual conditions
            trimTimeIdx = trialTimeIdx(1);
            % trim arrays of length = length(time)
            time = time(trimTimeIdx:end);
            Xdva = Xdva(trimTimeIdx:end);
            Ydva = Ydva(trimTimeIdx:end);
            if settings.driftCorrect
                dcX = dcX(trimTimeIdx:end);
                dcY = dcY(trimTimeIdx:end);
            end
            BSidx = BSidx(trimTimeIdx:end);
            fixIdx = fixIdx(trimTimeIdx:end);

            % set data arrays to max value of axes such that dot hovers at edge of
            % screen
            if settings.driftCorrect
                % X coordinates
                dcX(dcX>settings.axesLimit) = settings.axesLimit;
                dcX(dcX<-settings.axesLimit) = -settings.axesLimit;
                % Y coordinates
                dcY(dcY>settings.axesLimit) = settings.axesLimit;
                dcY(dcY<-settings.axesLimit) = -settings.axesLimit;
            else
                % X coordinates
                Xdva(Xdva>settings.axesLimit) = settings.axesLimit;
                Xdva(Xdva<-settings.axesLimit) = -settings.axesLimit;
                % Y coordinates
                Ydva(Ydva>settings.axesLimit) = settings.axesLimit;
                Ydva(Ydva<-settings.axesLimit) = -settings.axesLimit;
            end

            %% make movie

            % create figure
            figure
            hold off;

            %preallocate
            eventStr = 'no event';
            % M(1:length(time)+1) = struct('cdata',[],'colormap',[]);
            smallSaccCount = 0;
            bigSaccCount = 0;
            blinkCount = 0;
            fixCount = 0;
            movieFrame = 1;

            if settings.driftCorrect
                % set axis limits
                axis([-settings.axesLimit settings.axesLimit -settings.axesLimit settings.axesLimit]);

                % boundaries of screen for text placement
                xl = xlim; yl = ylim;
                xrange = xl(2)-xl(1); yrange = yl(2)-yl(1);

                % create movie frames
                M(movieFrame) = getframe;
                for iT = 1:length(time)

                    % check whether block has ended
                    if iT~=1 && any(trialStartTimes==time(iT))
                        % set inputs to AK_saveMovie()
                        clear part filename
                        part = ['Part' num2str(find(trialStartTimes==time(iT))-1)];
                        if settings.driftCorrect
                            filename = fullfile(saveParams.directory,[edfName '_DC_Run' num2str(runIdx) '_' part '.avi']);
                        else
                            filename = fullfile(saveParams.directory,[edfName '_noDC_Run' num2str(runIdx) '_' part '.avi']);
                        end
                        % save movie
                        AK_saveMovie(filename,M,saveParams.frameRate,saveParams.quality);

                        % create new movie (for memory reasons)
                        clear M
                        movieFrame = 1;
                        M(movieFrame) = getframe;
                    end

                    % advance movieFrame counter
                    movieFrame = movieFrame+1;

                    % choose color; encode blinks, saccades, fixations
                    clear C
                    if isnan(BSidx(iT))
                        C = [0 0 1];
                    elseif BSidx(iT) == 1
                        C = [1 0 0];
                    elseif BSidx(iT) == .75
                        C = [1 .8 0];
                    elseif fixIdx(iT) == 0 % (zero is fixation)
                        C = [0 .8 0];
                    else
                        C = [0 0 0];
                    end

                    % count small and large saccades and blinks
                    if iT~=1 && BSidx(iT)~=BSidx(iT-1) % has flagged classification of eye motion changed
                        if isnan(BSidx(iT)) && ~isnan(BSidx(iT-1)) % have to reiterate because nan~=nan is true
                            blinkCount = blinkCount+1;
                        elseif BSidx(iT) == 1
                            bigSaccCount = bigSaccCount+1;
                        elseif BSidx(iT) == .75^2 % <-- this is important: the BSidx for each saccade is multiplied by it's max(BSidx)
                            smallSaccCount = smallSaccCount+1;
                        end
                    end
                    % count n of times outside of fixation fixation
                    if iT~=1 && fixIdx(iT)~=fixIdx(iT-1) && fixIdx(iT)==1
                        fixCount = fixCount+1;
                    end

                    cla;

                    % stim zone
                    patch([-1 -1 1 1],[-1 1 1 -1],[1 1 1],'EdgeColor','k');
                    hold on;
                    % fixation zone
                    patch([-.5 -.5 .5 .5],[-.5 .5 .5 -.5],[1 1 1],'EdgeColor','k');

                    % plot xy-coordinates
                    scatter(dcX(iT),dcY(iT),'MarkerFaceColor',C,'MarkerEdgeColor','k'); % dc data
                    scatter(Xdva(iT),Ydva(iT),'MarkerFaceColor','w','MarkerEdgeColor','k'); % raw data
                    axis([-settings.axesLimit settings.axesLimit -settings.axesLimit settings.axesLimit]); % reiterate axis limits
                    xlabel('X (dva)'); ylabel('Y (dva)');
                    title(['Gaze Position Movie for ' edfName]);
                    hold off;

                    % time counter text box
                    text(xl(1)+50*xrange/1000,yl(1)+50*yrange/1000,['time = ' num2str((time(iT)-time(1))/1000)],... 
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');

                    % condition name text box
                    if any(time(iT) == trialStartTimes) % update event name
                        clear eventStr
                        eventStr = allEvents(trialStartTimes == time(iT));
                    end
                    text(xl(1)+50*xrange/1000,yl(2)-100*yrange/1000,['cond = ' eventStr],... 
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');

                    % blink and saccade text boxes
                    text(xl(2)-295*xrange/1000,yl(2)-50*yrange/1000,['blink count = ' num2str(blinkCount)],... % blinks
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');
                    text(xl(2)-420*xrange/1000,yl(2)-120*yrange/1000,['small saccade count = ' num2str(smallSaccCount)],... % small saccades
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');
                    text(xl(2)-390*xrange/1000,yl(2)-190*yrange/1000,['big saccade count = ' num2str(bigSaccCount)],... % big saccades
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');
                    % fixation text box
                    text(xl(2)-410*xrange/1000,yl(2)-260*yrange/1000,['fixation break count = ' num2str(fixCount)],... % breaks in fixation saccades
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');

                    M(movieFrame) = getframe;

                end
            else
                % set axis limits
                axis([-settings.axesLimit settings.axesLimit -settings.axesLimit settings.axesLimit]);

                % boundaries of screen for text placement
                xl = xlim; yl = ylim;
                xrange = xl(2)-xl(1); yrange = yl(2)-yl(1);

                % create movie frames
                M(movieFrame) = getframe;
                for iT = 1:length(time)

                    % check whether block has ended
                    if iT~=1 && any(trialStartTimes==time(iT))
                        % set inputs to AK_saveMovie()
                        clear part filename
                        part = ['Part' num2str(find(trialStartTimes==time(iT))-1)];
                        if settings.driftCorrect
                            filename = fullfile(saveParams.directory,[edfName '_DC_Run' num2str(runIdx) '_' part '.avi']);
                        else
                            filename = fullfile(saveParams.directory,[edfName '_noDC_Run' num2str(runIdx) '_' part '.avi']);
                        end
                        % save movie
                        AK_saveMovie(filename,M,saveParams.frameRate,saveParams.quality);

                        % create new movie (for memory reasons)
                        clear M
                        movieFrame = 1;
                        M(movieFrame) = getframe;
                    end

                    % advance movieFrame counter
                    movieFrame = movieFrame+1;

                    % choose color; encode blinks, saccades, fixations
                    clear C
                    if isnan(BSidx(iT))
                        C = [0 0 1];
                    elseif BSidx(iT) == 1
                        C = [1 0 0];
                    elseif BSidx(iT) == .75
                        C = [1 .8 0];
                    elseif fixIdx(iT) == 0 % (zero is fixation)
                        C = [0 .8 0];
                    else
                        C = [0 0 0];
                    end

                    % count small and large saccades and blinks
                    if iT~=1 && BSidx(iT)~=BSidx(iT-1) % has flagged classification of eye motion changed
                        if isnan(BSidx(iT)) && ~isnan(BSidx(iT-1)) % have to reiterate because nan~=nan is true
                            blinkCount = blinkCount+1;
                        elseif BSidx(iT) == 1
                            bigSaccCount = bigSaccCount+1;
                        elseif BSidx(iT) == .75^2 % <-- this is important: the BSidx for each saccade is multiplied by it's max(BSidx)
                            smallSaccCount = smallSaccCount+1;
                        end
                    end
                    % count n of times outside of fixation fixation
                    if iT~=1 && fixIdx(iT)~=fixIdx(iT-1) && fixIdx(iT)==1
                        fixCount = fixCount+1;
                    end

                    cla;

                    % stim zone
                    patch([-1 -1 1 1],[-1 1 1 -1],[1 1 1],'EdgeColor','k');
                    hold on;
                    % fixation zone
                    patch([-.5 -.5 .5 .5],[-.5 .5 .5 -.5],[1 1 1],'EdgeColor','k');

                    % plot xy-coordinates
                    scatter(Xdva(iT),Ydva(iT),'MarkerFaceColor',C,'MarkerEdgeColor','k');
                    axis([-settings.axesLimit settings.axesLimit -settings.axesLimit settings.axesLimit]); % reiterate axis limits
                    xlabel('X (dva)'); ylabel('Y (dva)');
                    title(['Gaze Position Movie for ' edfName]);
                    hold off;

                    % time counter text box
                    text(xl(1)+50*xrange/1000,yl(1)+50*yrange/1000,['time = ' num2str((time(iT)-time(1))/1000)],... 
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');

                    % condition name text box
                    if any(time(iT) == trialStartTimes) % update event name
                        clear eventStr
                        eventStr = allEvents(trialStartTimes == time(iT));
                    end
                    text(xl(1)+50*xrange/1000,yl(2)-100*yrange/1000,['cond = ' eventStr],... 
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');

                    % blink and saccade text boxes
                    text(xl(2)-295*xrange/1000,yl(2)-50*yrange/1000,['blink count = ' num2str(blinkCount)],... % blinks
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');
                    text(xl(2)-420*xrange/1000,yl(2)-120*yrange/1000,['small saccade count = ' num2str(smallSaccCount)],... % small saccades
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');
                    text(xl(2)-390*xrange/1000,yl(2)-190*yrange/1000,['big saccade count = ' num2str(bigSaccCount)],... % big saccades
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');
                    % fixation text box
                    text(xl(2)-410*xrange/1000,yl(2)-260*yrange/1000,['fixation break count = ' num2str(fixCount)],... % breaks in fixation saccades
                        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                        'hittest','off');

                    M(movieFrame) = getframe;

                end
            end

            % save last movie part:
            % set inputs to AK_saveMovie()
            clear part filename
            part = ['Part' num2str(length(trialStartTimes))];
            if settings.driftCorrect
                filename = fullfile(saveParams.directory,[edfName '_DC_Run' num2str(runIdx) '_' part '.avi']);
            else
                filename = fullfile(saveParams.directory,[edfName '_noDC_Run' num2str(runIdx) '_' part '.avi']);
            end
            % save movie
            AK_saveMovie(filename,M,saveParams.frameRate,saveParams.quality);
        end
    end
else
    % load raw data
    if ~exist(mat_filename,'file') % check whether or not function output has already been saved
        disp(['parsing ' asc_filename]) % message
        [timepts, event] = AK_GABAinASD_ascfileParse(asc_filename); % parse asc file
        save(mat_filename,'timepts','event'); % save function output
    else
        disp(['loading ' mat_filename]) % message
        load(mat_filename,'timepts','event') % load file containing timepts array
    end
    
    %% processing: create indices flagging fixation, saccades, and blinks; create indices marking conditions and extract list of condition names; drift correct

    % extract and prepare raw data from timepts
    time = cellfun(@str2double,timepts(2:end,1)); % set time vector
    Xpixel = cellfun(@str2double,timepts(2:end,2)); % set vector of X position in pixels
    Ypixel = cellfun(@str2double,timepts(2:end,3)); % set vector of Y position in pixels
    Xdva = (Xpixel-display.center(1))./pixelsPerDegree(1); % set vector of X position as degrees visual angle
    Ydva = (Ypixel-display.center(2))./pixelsPerDegree(2); % set vector of Y position as degrees visual angle
    events = timepts(2:length(timepts),6); % store list of events as cell array of strings

    % blinks and saccades:
    % find number of blinks and idx into time array where there are blinks
    blink.N = length(event.blink(:,1))-1;
    isBlink = zeros(length(time),1);
    if blink.N>0 % is there at least one blink?
        for bl = 2:length(event.blink(:,1)); % cycle through saccades
            % get blink start times and end times
            blink.Stimes(bl-1) = str2double(event.blink{bl,1});
            blink.Etimes(bl-1) = str2double(event.blink{bl,2});
            % find an index in time for each blink
            isBlink(blink.Stimes(bl-1)<= time & time<=blink.Etimes(bl-1)) = 1;
        end
    else % if there are no blinks
        blink.Stimes = nan;
        blink.Etimes = nan;
    end
    % find number of saccades and idx into time array where there are saccades
    sacc.N = length(event.sacc(:,1))-1;
    isSacc = zeros(length(time),1);
    if sacc.N>0 % is there at least one saccade?
        for sa = 2:length(event.sacc(:,1)); % cycle through saccades
            % get saccade start times and end times
            sacc.Stimes(sa-1) = str2double(event.sacc{sa,1});
            sacc.Etimes(sa-1) = str2double(event.sacc{sa,2});
            sacc.Amps(sa-1) = str2double(event.sacc{sa,8});
            % find an index in time for each saccade
            if sacc.Amps(sa-1) > settings.fixationRadiusThreshold % for big saccades
                isSacc(sacc.Stimes(sa-1)<= time & time<=sacc.Etimes(sa-1)) = 1;
            else % for small saccades
                isSacc(sacc.Stimes(sa-1)<= time & time<=sacc.Etimes(sa-1)) = .75;
            end
        end
    else % if there are no saccades
        sacc.Stimes = nan;
        sacc.Etimes = nan;
    end
    % differentiate small saccades, big saccades and blinks
    BSidx = isBlink+isSacc;
    if length(find(BSidx==0)) > 1
        for iSa = 1:sacc.N
            clear tempIdxS tempIdxE
            tempIdxS = find(time==sacc.Stimes(iSa)); tempIdxE = find(time==sacc.Etimes(iSa));
            BSidx(tempIdxS:tempIdxE) = BSidx(tempIdxS:tempIdxE)*max(BSidx(tempIdxS:tempIdxE));
        end
        BSidx(BSidx>1) = nan; % use to distinguish blinks from saccades
    end

    % segment information about events and create indices
    clear eventChangeIdx allEvents eventTimes
    eventChangeIdx = find(strcmp(events(1:end-1),events(2:end))==0); % create index of when events change
    allEvents = events(eventChangeIdx); % create list of events in order
    eventTimes = time(eventChangeIdx); % create vector of time points at which events change

    % drift correction of coordinates and fixation index
    % preallocate
    fixIdx = zeros(length(time),1);
    if settings.driftCorrect
        % NEW DRIFT CORRECTION: linear transformations optimized for each block based on fixation coordinates
        % define what counts as a fixation
        fixations = AK_defineFixations([Xdva,Ydva,time],settings.fixationRadiusThreshold,settings.fixationMinDur);

        % separate fixations into chunks by block and drift correct within each chunk
        dcX = []; % preallocate variables
        dcY = [];
        chunkBoundIdx = [1;eventChangeIdx;length(time)]; % beginning and end indices in time of each chunk
        for iChunk = 1:length(chunkBoundIdx)-1
            useFix = zeros(length(fixations),1); % reset index of fixations to use for this chunk
            % determine which fixations to use in the current chunk's drift correction
            for iFix = 1:length(fixations)
                if ~isempty(intersect(fixations(iFix).timeIdx,chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1))) % compare indices for fixation to indices for chunk
                    useFix(iFix) = 1; % use this fixation for drift correction in the current chunk
                end    
            end
            useFix = logical(useFix);
            % prepare drift correction function input
            clear chunkFix chunkData chunkDCcoords
            chunkFix = vertcat(fixations(useFix).coordinates)'; % fixation coordinates for this chunk (n fixations by x & y)
            chunkData = horzcat(Xdva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)),Ydva(chunkBoundIdx(iChunk):chunkBoundIdx(iChunk+1)))';

            % apply drift correction to chunk
            chunkDCcoords = AK_driftCorrect(chunkData,chunkFix,settings.stimCoords);
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

        % create fixation idx (zero is fixation)
        fixIdx(sqrt(dcX.^2 + dcY.^2) > settings.fixationRadiusThreshold) = 1;
    else
        % create fixation idx (zero is fixation)
        fixIdx(sqrt(Xdva.^2 + Ydva.^2) > settings.fixationRadiusThreshold) = 1;
    end

    % trim data prior to actual conditions
    trimTimeIdx = find(time==eventTimes(3));
    % trim arrays of length = length(time)
    time = time(trimTimeIdx:end);
    Xdva = Xdva(trimTimeIdx:end);
    Ydva = Ydva(trimTimeIdx:end);
    if settings.driftCorrect
        dcX = dcX(trimTimeIdx:end);
        dcY = dcY(trimTimeIdx:end);
    end
    BSidx = BSidx(trimTimeIdx:end);
    fixIdx = fixIdx(trimTimeIdx:end);
    % trim event arrays
    allEvents = allEvents(3:end);
    eventTimes = eventTimes(3:end);

    % set data arrays to max value of axes such that dot hovers at edge of
    % screen
    if settings.driftCorrect
        % X coordinates
        dcX(dcX>settings.axesLimit) = settings.axesLimit;
        dcX(dcX<-settings.axesLimit) = -settings.axesLimit;
        % Y coordinates
        dcY(dcY>settings.axesLimit) = settings.axesLimit;
        dcY(dcY<-settings.axesLimit) = -settings.axesLimit;
    else
        % X coordinates
        Xdva(Xdva>settings.axesLimit) = settings.axesLimit;
        Xdva(Xdva<-settings.axesLimit) = -settings.axesLimit;
        % Y coordinates
        Ydva(Ydva>settings.axesLimit) = settings.axesLimit;
        Ydva(Ydva<-settings.axesLimit) = -settings.axesLimit;
    end

    %% make movie

    % create figure
    figure
    hold off;

    %preallocate
    eventStr = 'no event';
    % M(1:length(time)+1) = struct('cdata',[],'colormap',[]);
    smallSaccCount = 0;
    bigSaccCount = 0;
    blinkCount = 0;
    fixCount = 0;
    movieFrame = 1;

    if settings.driftCorrect
        % set axis limits
        axis([-settings.axesLimit settings.axesLimit -settings.axesLimit settings.axesLimit]);

        % boundaries of screen for text placement
        xl = xlim; yl = ylim;
        xrange = xl(2)-xl(1); yrange = yl(2)-yl(1);

        % create movie frames
        M(movieFrame) = getframe;
        for iT = 1:length(time)

            % check whether block has ended
            if iT~=1 && any(eventTimes==time(iT))
                % set inputs to AK_saveMovie()
                clear part filename
                part = ['Part' num2str(find(eventTimes==time(iT))-1)];
                if settings.driftCorrect
                    filename = fullfile(saveParams.directory,[edfName '_DC_' part '.avi']);
                else
                    filename = fullfile(saveParams.directory,[edfName '_noDC_' part '.avi']);
                end
                % save movie
                AK_saveMovie(filename,M,saveParams.frameRate,saveParams.quality);

                % create new movie (for memory reasons)
                clear M
                movieFrame = 1;
                M(movieFrame) = getframe;
            end

            % advance movieFrame counter
            movieFrame = movieFrame+1;

            % choose color; encode blinks, saccades, fixations
            clear C
            if isnan(BSidx(iT))
                C = [0 0 1];
            elseif BSidx(iT) == 1
                C = [1 0 0];
            elseif BSidx(iT) == .75
                C = [1 .8 0];
            elseif fixIdx(iT) == 0 % (zero is fixation)
                C = [0 .8 0];
            else
                C = [0 0 0];
            end

            % count small and large saccades and blinks
            if iT~=1 && BSidx(iT)~=BSidx(iT-1) % has flagged classification of eye motion changed
                if isnan(BSidx(iT)) && ~isnan(BSidx(iT-1)) % have to reiterate because nan~=nan is true
                    blinkCount = blinkCount+1;
                elseif BSidx(iT) == 1
                    bigSaccCount = bigSaccCount+1;
                elseif BSidx(iT) == .75^2 % <-- this is important: the BSidx for each saccade is multiplied by it's max(BSidx)
                    smallSaccCount = smallSaccCount+1;
                end
            end
            % count n of times outside of fixation fixation
            if iT~=1 && fixIdx(iT)~=fixIdx(iT-1) && fixIdx(iT)==1
                fixCount = fixCount+1;
            end

            cla;

            % stim zone
            patch([-1 -1 1 1],[-1 1 1 -1],[1 1 1],'EdgeColor','k');
            hold on;
            % fixation zone
            patch([-.5 -.5 .5 .5],[-.5 .5 .5 -.5],[1 1 1],'EdgeColor','k');

            % plot xy-coordinates
            scatter(dcX(iT),dcY(iT),'MarkerFaceColor',C,'MarkerEdgeColor','k'); % dc data
            scatter(Xdva(iT),Ydva(iT),'MarkerFaceColor','w','MarkerEdgeColor','k'); % raw data
            axis([-settings.axesLimit settings.axesLimit -settings.axesLimit settings.axesLimit]); % reiterate axis limits
            xlabel('X (dva)'); ylabel('Y (dva)');
            title(['Gaze Position Movie for ' edfName]);
            hold off;

            % time counter text box
            text(xl(1)+50*xrange/1000,yl(1)+50*yrange/1000,['time = ' num2str((time(iT)-time(1))/1000)],... 
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');

            % condition name text box
            if any(time(iT) == eventTimes) % update event name
                clear eventStr
                eventStr = allEvents(eventTimes == time(iT));
            end
            text(xl(1)+50*xrange/1000,yl(2)-100*yrange/1000,['cond = ' eventStr],... 
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');

            % blink and saccade text boxes
            text(xl(2)-295*xrange/1000,yl(2)-50*yrange/1000,['blink count = ' num2str(blinkCount)],... % blinks
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');
            text(xl(2)-420*xrange/1000,yl(2)-120*yrange/1000,['small saccade count = ' num2str(smallSaccCount)],... % small saccades
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');
            text(xl(2)-390*xrange/1000,yl(2)-190*yrange/1000,['big saccade count = ' num2str(bigSaccCount)],... % big saccades
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');
            % fixation text box
            text(xl(2)-410*xrange/1000,yl(2)-260*yrange/1000,['fixation break count = ' num2str(fixCount)],... % breaks in fixation saccades
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');

            M(movieFrame) = getframe;

        end
    else
        % set axis limits
        axis([-settings.axesLimit settings.axesLimit -settings.axesLimit settings.axesLimit]);

        % boundaries of screen for text placement
        xl = xlim; yl = ylim;
        xrange = xl(2)-xl(1); yrange = yl(2)-yl(1);

        % create movie frames
        M(movieFrame) = getframe;
        for iT = 1:length(time)

            % check whether block has ended
            if iT~=1 && any(eventTimes==time(iT))
                % set inputs to AK_saveMovie()
                clear part filename
                part = ['Part' num2str(find(eventTimes==time(iT))-1)];
                if settings.driftCorrect
                    filename = fullfile(saveParams.directory,[edfName '_DC_' part '.avi']);
                else
                    filename = fullfile(saveParams.directory,[edfName '_noDC_' part '.avi']);
                end
                % save movie
                AK_saveMovie(filename,M,saveParams.frameRate,saveParams.quality);

                % create new movie (for memory reasons)
                clear M
                movieFrame = 1;
                M(movieFrame) = getframe;
            end

            % advance movieFrame counter
            movieFrame = movieFrame+1;

            % choose color; encode blinks, saccades, fixations
            clear C
            if isnan(BSidx(iT))
                C = [0 0 1];
            elseif BSidx(iT) == 1
                C = [1 0 0];
            elseif BSidx(iT) == .75
                C = [1 .8 0];
            elseif fixIdx(iT) == 0 % (zero is fixation)
                C = [0 .8 0];
            else
                C = [0 0 0];
            end

            % count small and large saccades and blinks
            if iT~=1 && BSidx(iT)~=BSidx(iT-1) % has flagged classification of eye motion changed
                if isnan(BSidx(iT)) && ~isnan(BSidx(iT-1)) % have to reiterate because nan~=nan is true
                    blinkCount = blinkCount+1;
                elseif BSidx(iT) == 1
                    bigSaccCount = bigSaccCount+1;
                elseif BSidx(iT) == .75^2 % <-- this is important: the BSidx for each saccade is multiplied by it's max(BSidx)
                    smallSaccCount = smallSaccCount+1;
                end
            end
            % count n of times outside of fixation fixation
            if iT~=1 && fixIdx(iT)~=fixIdx(iT-1) && fixIdx(iT)==1
                fixCount = fixCount+1;
            end

            cla;

            % stim zone
            patch([-1 -1 1 1],[-1 1 1 -1],[1 1 1],'EdgeColor','k');
            hold on;
            % fixation zone
            patch([-.5 -.5 .5 .5],[-.5 .5 .5 -.5],[1 1 1],'EdgeColor','k');

            % plot xy-coordinates
            scatter(Xdva(iT),Ydva(iT),'MarkerFaceColor',C,'MarkerEdgeColor','k');
            axis([-settings.axesLimit settings.axesLimit -settings.axesLimit settings.axesLimit]); % reiterate axis limits
            xlabel('X (dva)'); ylabel('Y (dva)');
            title(['Gaze Position Movie for ' edfName]);
            hold off;

            % time counter text box
            text(xl(1)+50*xrange/1000,yl(1)+50*yrange/1000,['time = ' num2str((time(iT)-time(1))/1000)],... 
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');

            % condition name text box
            if any(time(iT) == eventTimes) % update event name
                clear eventStr
                eventStr = allEvents(eventTimes == time(iT));
            end
            text(xl(1)+50*xrange/1000,yl(2)-100*yrange/1000,['cond = ' eventStr],... 
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');

            % blink and saccade text boxes
            text(xl(2)-295*xrange/1000,yl(2)-50*yrange/1000,['blink count = ' num2str(blinkCount)],... % blinks
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');
            text(xl(2)-420*xrange/1000,yl(2)-120*yrange/1000,['small saccade count = ' num2str(smallSaccCount)],... % small saccades
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');
            text(xl(2)-390*xrange/1000,yl(2)-190*yrange/1000,['big saccade count = ' num2str(bigSaccCount)],... % big saccades
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');
            % fixation text box
            text(xl(2)-410*xrange/1000,yl(2)-260*yrange/1000,['fixation break count = ' num2str(fixCount)],... % breaks in fixation saccades
                'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                'hittest','off');

            M(movieFrame) = getframe;

        end
    end

    % save last movie part:
    % set inputs to AK_saveMovie()
    clear part filename
    part = ['Part' num2str(length(eventTimes))];
    if settings.driftCorrect
        filename = fullfile(saveParams.directory,[edfName '_DC_' part '.avi']);
    else
        filename = fullfile(saveParams.directory,[edfName '_noDC_' part '.avi']);
    end
    % save movie
    AK_saveMovie(filename,M,saveParams.frameRate,saveParams.quality);
end

end

