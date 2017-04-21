function [ output_args ] = AK_GABA_getPupilSizeData( ascFilename, directory )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% check input
if nargin < 1
    error('AK_GABA_getPupilSizeData requires at least an asc filename (string) as an argument');
end
if nargin < 2
    directory = cd;
end

%% construct directory and load data

mat_name = [fullfile(directory,ascFilename(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
if ~exist(mat_name,'file') % check whether or not function output has already been saved
    disp(['parsing ' fullfile(directory,ascFilename)]) % message
    block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(directory,ascFilename)); % parse asc file
    save(mat_name,'block'); % save function output
else
    disp(['loading ' mat_name]) % message
    load(mat_name,'block') % load file containing timepts array
end
%% determine whether or not to correct for eye link message timing bug from early version of the experiment based on the date each edf was written

correctionDate = datetime('2015-12-17','InputFormat','yyyy-MM-dd'); % the day after the last participant was run on the uncorrected code
avgError = 997; % the average error in the timing of the eyetracker messages for stimulus events

% find analogous .edf for this .asc file
clear thisEdf correction hIdx
thisEdf = dir(fullfile(directory,[asc_list(iA).name(1:end-3) 'edf']));
% determine whether correction is needed by checking .edf date
correction = correctionDate > datetime(thisEdf.date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
% test: correction
disp(['analogous .edf name = ' thisEdf.name '; file date = ' thisEdf.date]);
disp(['apply correction? ' num2str(correction)]);
% end test: correction

%% process pupil size data by extracting raw signal and baseline measures

for iB = 1:length(block); % cycle through blocks of trials
    % only proceed if there is a full set of events and time is continuous    
    if length(block(iB).events) > 21 && ~any((cellfun(@str2double,block(iB).timepts(3:end,1))-cellfun(@str2double,block(iB).timepts(2:end-1,1)))~=1) 
        clear time pSize pSizeFilt events eventTimes blink Stim1Index Stim2Index timeIndexS  % clear variables

        % extract relevant data from block struct
        time{iS,iB} = cellfun(@str2double,block(iB).timepts(2:end,1)); % set time vector
        pSize = cellfun(@str2double,block(iB).timepts(2:end,4)); % set pupil size vector

        events = block(iB).events(2:end,2); % store list of events as cell array of strings
        eventTimes = cellfun(@str2double,block(iB).events(2:end,1)); % set vector of event times
        events(cellfun(@isempty,events)) = {'no event'}; % fill in empty events

        if removeBlinks
            % get blink times
            blink = [];
            for iBl = 2:length(block(iB).blink(:,1))
                blink = [blink str2double(block(iB).blink{iBl,1}):str2double(block(iB).blink{iBl,2})];
            end 
            % get saccade times
            sacc = []; 
            for iSa = 2:length(block(iB).sacc(:,1))
                sacc = [sacc str2double(block(iB).sacc{iSa,1}):str2double(block(iB).sacc{iSa,2})];
            end 
        end

        % create indices for event transitions
        Stim1Index = cellfun(@regexp,events,repmat({regexptranslate('wildcard','*Stim1*')},size(events)),'UniformOutput',0); % index of stim1 events
        Stim1Index(cellfun(@isempty,Stim1Index)) = {0}; % fill in empty cells
        Stim1Index = logical(cell2mat(Stim1Index)); % Convert to logical for indexing
        Stim2Index = cellfun(@regexp,events,repmat({regexptranslate('wildcard','*Stim2*')},size(events)),'UniformOutput',0); % index of stim2 events
        Stim2Index(cellfun(@isempty,Stim2Index)) = {0}; % fill in empty cells
        Stim2Index = logical(cell2mat(Stim2Index)); % Convert to logical for indexing

        if correction == 1 % correct for error in stimulus time coding
            WhiteOn{iS,iB} = eventTimes(Stim2Index)+avgError; % time of interest start times: black to white
            BlackOn{iS,iB} = eventTimes(Stim1Index)+avgError; % white to black
        else
            WhiteOn{iS,iB} = eventTimes(Stim2Index); % time of interest start times: black to white
            BlackOn{iS,iB} = eventTimes(Stim1Index); % white to black
        end

        % remove far outliers from pupil data using IQR filter
        switch filter
            case 'IQR'
                % AMK IQR outliers->nan filter
                [pSize,oIdx{iB}] = AK_windowIQRfilt(pSize,4000,3); % eliminate outliers outside of the IQR by +/- 1.5*IQR, w/in a sliding 2000ms window  
            case 'Hampel'
                % Hampel outliers->windowMedian filter
                [pSize,oIdx{iB}] = hampel(pSize,2000,3);
            case 'SG'
                % apply Savitzky-Golay smoothing filter; local cubic modeling w/in 51 frame window 
                pSize = sgolayfilt(pSize,3,51); 
                  % modified frame rate from Maria's code, where param = 13:
                  % 51=13*(1000hz/250hz)-1 (frame rate must be odd)
            case 'windowAvg'
                % AMK window average smoothing filter
                pSize = AK_windowAvgFilt(pSize,400,avgFunction);
            case 'BLD'
                % AMK band-limited differentiator filter (Low-Pass)
                [~,pSize,~] = AK_BLDfilt(pSize,.5,.6,1000); % bandpass = .5; bandstop = 1; hz = 1000
        end

        % subtract out baseline
        baseline{iS,iA}(iB) = avgFunction(pSize);
        pSizeCell{iB} = pSize' - baseline{iS,iA}(iB);

        % trim pupil data timecourse to start at the first
        % BlackOn and end at last BlackOn
%                         pSizeCell{iB}(find(time{iS,iB}==BlackOn{iS,iB}(end)):end) = []; % trim end first to avoid indexing error
        pSizeCell{iB}(1:find(time{iS,iB}==BlackOn{iS,iB}(1))-1) = [];
        if any(strcmp(filter,{'IQR','Hampel'}))
            % trim index of points affected by filtering
%                             hIdx{iB}(find(time{iS,iB}==BlackOn{iS,iB}(end)):end) = []; 
            oIdx{iB}(1:find(time{iS,iB}==BlackOn{iS,iB}(1))-1) = [];
        end
%                         time{iS,iB}(find(time{iS,iB}==BlackOn{iS,iB}(end)):end) = []; % trim time last so that everything is trimmed the same
        time{iS,iB}(1:find(time{iS,iB}==BlackOn{iS,iB}(1))-1) = [];

        % index in time where events start
        timeIndexS = arrayfun(@(x) find(time{iS,iB}==x),WhiteOn{iS,iB});

        for TI = 1:length(WhiteOn{iS,iB}); % cycle through times of interest: times when the white stimulus comes on  
            if removeBlinks
                % convert blink and saccade trials to nans
                if TI~=length(WhiteOn{iS,iB}) && (~isempty(intersect(blink,(WhiteOn{iS,iB}(TI)-500:WhiteOn{iS,iB}(TI+1)))) || ~isempty(intersect(sacc,(WhiteOn{iS,iB}(TI)-500:WhiteOn{iS,iB}(TI+1))))) 
                    pSizeCell{iB}(timeIndexS(TI):timeIndexS(TI+1)-1) = nan;
                elseif ~isempty(intersect(blink,(WhiteOn{iS,iB}(TI)-500:time{iS,iB}(end)))) || ~isempty(intersect(sacc,(WhiteOn{iS,iB}(TI)-500:time{iS,iB}(end))))
                    pSizeCell{iB}(timeIndexS(TI):end) = nan;
                end
            end
            if strcmp(plotWhat,'GroupAvgTrace')
                % relative pupil size trace for each trial
                try
                    pSizeTrials{iS}(TI+10*(iB-1),:) = pSizeCell{iB}(timeIndexS(TI):timeIndexS(TI)+1999) - pSizeCell{iB}(timeIndexS(TI)-250);
                catch
                    if ~isempty(pSizeCell{iB}(timeIndexS(TI):end))
                        pSizeTrials{iS}(TI+10*(iB-1),:) = [pSizeCell{iB}(timeIndexS(TI):end) nan(1,2000-length(pSizeCell{iB}(timeIndexS(TI):end)))] - pSizeCell{iB}(timeIndexS(TI)-250);
                    else
                        pSizeTrials{iS}(TI+10*(iB-1),:) = nan(1,2000);
                    end
                end
            end
        end

%                         % wtf is going on?
%                         test(iS,iA,iB).length = length(pSize);
%                         test(iS,iA,iB).nWhite = length(WhiteOn{iS,iB});
%                         test(iS,iA,iB).nBlack = length(BlackOn{iS,iB});
%                         test(iS,iA,iB).eventNames = events;
%                         test(iS,iA,iB).eventTimes = eventTimes;
%                         test(iS,iA,iB).whiteIdx = Stim2Index;
%                         test(iS,iA,iB).blackIdx = Stim1Index;
    end
end

% set up figures for traces per subject
if strcmp(plotWhat,'trace')
    figure(iS);hold on;
    blockColors = {'r','g','m','c'};
end

% nan pad to make each block the same length for concatenation
%                 maxBlockLength = max(cellfun(@length,pSizeCell));
maxBlockLength = 22000;
for iB = 1:length(pSizeCell)
    if any(size(pSizeCell{iB})==0) % where there is no data, everything is nans
        pSizeCell{iB} = nan([1 maxBlockLength]);
        if any(strcmp(filter,{'IQR','Hampel'}))
            oIdx{iB} = zeros([1 maxBlockLength]);
        end
    else
        pSizeCell{iB} = AK_nanPad(pSizeCell{iB},[1 maxBlockLength],[1 1]);
        if any(strcmp(filter,{'IQR','Hampel'}))
            oIdx{iB} = AK_nanPad(oIdx{iB},[1 maxBlockLength],[1 1]);
            oIdx{iB}(isnan(oIdx{iB})) = 0;
        end
    end
    % trim each array to equal length
    pSizeCell{iB} = pSizeCell{iB}(1:maxBlockLength);
    if any(strcmp(filter,{'IQR','Hampel'}))
        oIdx{iB} = oIdx{iB}(1:maxBlockLength);
    end
    if strcmp(plotWhat,'trace')
        % plot trace for each block along with avg (SG filter to smoothe)
        if ~strcmp(filter,'SG') || (~strcmp(filter,'windowAvg') && ~isequal(avgFunction,@nanmean))
            pSizeFilt{iB} = sgolayfilt(pSizeCell{iB},3,51);
            plot(1:length(pSizeFilt{iB}),pSizeFilt{iB},['--' blockColors{iB}]);
        else
            plot(1:length(pSizeCell{iB}),pSizeCell{iB},['--' blockColors{iB}]);
        end
        % plot event code loci
        for iW = 1:length(WhiteOn{iS,iB})
            plot([WhiteOn{iS,iB}(iW)-time{iS,iB}(1) WhiteOn{iS,iB}(iW)-time{iS,iB}(1)],[min(cell2mat(cellfun(@min,pSizeCell,'UniformOutput',false))) max(cell2mat(cellfun(@max,pSizeCell,'UniformOutput',false)))],'-.k');
        end
        for iBl = 1:length(BlackOn{iS,iB})
            plot([BlackOn{iS,iB}(iBl)-time{iS,iB}(1) BlackOn{iS,iB}(iBl)-time{iS,iB}(1)],[min(cell2mat(cellfun(@min,pSizeCell,'UniformOutput',false))) max(cell2mat(cellfun(@max,pSizeCell,'UniformOutput',false)))],':k');
        end
    end
end

% cat
clear pupilData 
pupilData = vertcat(pSizeCell{:});

% save data matrix for Scott
save(fullfile(home_dir,[subj_dirs{iS} '_PupilSizeRaw.mat']),'pupilData');
save(fullfile(home_dir,[subj_dirs{iS} '_PupilSizeFilt.mat']),'pupilData');


end

