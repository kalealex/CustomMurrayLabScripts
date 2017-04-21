% makeEyetrackingMovie_GABAinASD
% MurrayLab_2016
% Created by AMK 8/10/16
% new drift correction added 8/24/16


%% input

% user selections:
% file selection
subject = 'G114'; 
fMRIsession = 'fMRI2';
edfName = 'G114_16';
% drift correction params
driftCorrect = 1; % logical (do drift correction?)
fixationRadiusThreshold = 1; % used to be .5
fixationMinDur = 100; % fixations must be at least this many milliseconds
stimCoords = [0;0]; % coordinates of stimuli on screen (x & y by number of stim)
% axis limits
AxLim = 5; % display boundaries as distance from fixation (dva)
% movie params
save_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab\Eye_Tracking_Data\G114_EyetrackingMovie';
frameRate = 60; % defaults to 30 hz
quality = 50; % defaults to 75

%% enter display characteristics for calculating visual angle, pixels per degree

display.dist = 66; % distance to screen in cm
display.width = 33; % width of display in cm
display.height = 24.75; % height of display in cm
display.resolution = [1024 768]; % number of pixels for screen width, height
display.center = display.resolution./2; % screen center in pixel coordinates
display.angles = [2*atand((display.width/2)/display.dist) 2*atand((display.height/2)/display.dist)]; % visual angle of the screen width, height
pixelsPerDegree = display.resolution./display.angles; % pixels per degree visual angle for X dimension and Y dimension

%% make directory and load file

% create directory and file name
top_dir = 'L:\MurrayLab\ASD\Data';
data_dir = fullfile(top_dir,subject,fMRIsession,'log_files');
mat_filename = fullfile(data_dir,[edfName '.mat']);
asc_filename = fullfile(data_dir,[edfName '.asc']);

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
        if sacc.Amps(sa-1) > fixationRadiusThreshold % for big saccades
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
clear patchIndex patchEvents uniqueEvents patchTimes patchX patchY patchC
eventChangeIdx = find(strcmp(events(1:end-1),events(2:end))==0); % create index of when events change
allEvents = events(eventChangeIdx); % create list of events in order
eventTimes = time(eventChangeIdx); % create vector of time points at which events change

% drift correction of coordinates and fixation index
% preallocate
fixIdx = zeros(length(time),1);
if driftCorrect
    % NEW DRIFT CORRECTION: linear transformations optimized for each block based on fixation coordinates
    % define what counts as a fixation
    fixations = AK_defineFixations([Xdva,Ydva,time],fixationRadiusThreshold,fixationMinDur);
    
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
    
    % create fixation idx (zero is fixation)
    fixIdx(sqrt(dcX.^2 + dcY.^2) > fixationRadiusThreshold) = 1;
    
%     % OLD DRIFT CORRECTION: linear model to regress out drift; exclude blinks
%     try
%         % run linear model for X and Y dva
%         Xmdl = LinearModel.fit(time,Xdva,'linear','Exclude',isnan(Xdva));
%         Ymdl = LinearModel.fit(time,Ydva,'linear','Exclude',isnan(Ydva));
%         % store vector of residuals
%         dcX = Xmdl.Residuals.Raw;
%         dcY = Ymdl.Residuals.Raw;
%         % create fixation idx (zero is fixation)
%         fixIdx(sqrt(dcX.^2 + dcY.^2) > fixationThreshold) = 1;
%     catch
%         dcX = nan(length(Xdva),1);
%         dcY = nan(length(Ydva),1);
%     end
else
    % create fixation idx (zero is fixation)
    fixIdx(sqrt(Xdva.^2 + Ydva.^2) > fixationRadiusThreshold) = 1;
end

% trim data prior to actual conditions
trimTimeIdx = find(time==eventTimes(3));
% trim arrays of length = length(time)
time = time(trimTimeIdx:end);
Xdva = Xdva(trimTimeIdx:end);
Ydva = Ydva(trimTimeIdx:end);
if driftCorrect
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
if driftCorrect
    % X coordinates
    dcX(dcX>AxLim) = AxLim;
    dcX(dcX<-AxLim) = -AxLim;
    % Y coordinates
    dcY(dcY>AxLim) = AxLim;
    dcY(dcY<-AxLim) = -AxLim;
else
    % X coordinates
    Xdva(Xdva>AxLim) = AxLim;
    Xdva(Xdva<-AxLim) = -AxLim;
    % Y coordinates
    Ydva(Ydva>AxLim) = AxLim;
    Ydva(Ydva<-AxLim) = -AxLim;
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

switch driftCorrect
    case 1
        % set axis limits
%         minX = min(dcX); minY = min(dcY);
%         maxX = max(dcX); maxY = max(dcY);
%         AxLim = max([abs(minX) abs(minY) maxX maxY]);
        axis([-AxLim AxLim -AxLim AxLim]);
        
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
                filename = fullfile(save_dir,[edfName 'newDC2_' part '.avi']);
                
                % save movie
                AK_saveMovie(filename,M,frameRate,quality);
                
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
            axis([-AxLim AxLim -AxLim AxLim]); % reiterate axis limits
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
%         imagesc(gamma(128*ones(size(gratings{iSize}(:,:,1)))),[0 255]); 
%         M(iT+2) = getframe;
    case 0
        % set axis limits
%         minX = min(Xdva); minY = min(Ydva);
%         maxX = max(Xdva); maxY = max(Ydva);
%         AxLim = max([abs(minX) abs(minY) maxX maxY]);
        axis([-AxLim AxLim -AxLim AxLim]);
        
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
                filename = fullfile(save_dir,[edfName 'newDC2_' part '.avi']);
                
                % save movie
                AK_saveMovie(filename,M,frameRate,quality);
                
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
            axis([-AxLim AxLim -AxLim AxLim]); % reiterate axis limits
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
%         imagesc(gamma(128*ones(size(gratings{iSize}(:,:,1)))),[0 255]); 
%         M(iT+2) = getframe;
end

% save last movie part:
% set inputs to AK_saveMovie()
clear part filename
part = ['Part' num2str(length(eventTimes))];
filename = fullfile(save_dir,[edfName 'newDC2_' part '.avi']);

% save movie
AK_saveMovie(filename,M,frameRate,quality);


% % save movie w/ settings
% myMovie = VideoWriter(fullfile(save_dir,[edfName '_' part '.avi']));
% myMovie.FrameRate = 60; % Default 30
% myMovie.Quality = 50; % Default 75
% % write video
% open(myMovie);
% writeVideo(myMovie,M);
% close(myMovie);

% movie2avi(M,fullfile(data_dir,[edfName '.avi']),'fps',60)


