function [ trialPupilSize, trialBaseline, blockPupilSize, blockBaseline, avgPupilSizeAcrossBlocks, power ] = AK_GABA_getPupilSizeData( ascFilename,...
        directory, removeBlinks, avgFunction, filter, plotWhat )
%AK_GABA_getPupilSizeData retrieves and formats data from the GABAinASD
%Pupil Size experiment, in which a participant fixates while a disk
%oscillates from black to white with a period of two seconds. The full
%experiment is four blocks of 10 cycles each, for a total of 40 trials.
%This function gets the pupil size timecourse data and groups it by trial,
%by block, and as an overall average timecourse. The function also performs
%a Fourier power analysis. Additionally, it offers the options to remove
%blinks from the data, select a function for averaging, select from a
%variety of filters, and select from a variety of ways of visualizing the
%data.
%   INPUT:
%       ascFilename: an .asc file name as a string (i.e., 'G101p1.asc')
%       directory[optional]: the directory where the .asc lives (defaults to cd)
%       removeBlinks[optional]: a boolean (true => replace data during
%           blinks with nans; false => leave blinks in the data); defaults to
%           false
%       avgFunction[optional]: the handle of the function that is used for
%           all averaging during data processing (this includes the
%           calculation of block-wise baselines- which are subtracted from
%           the pupil timecourse of each block, the averaging of data
%           within and across blocks, and the averaging that is used in the
%           'windowAvg' filter); defaults to @nanmedian; either @nanmean or
%           @nanmedian are recommended since they filter out nans (blinks)
%       filter[optional]: choose a filter to apply to each block of data
%           from among the following options (copy the string exactly):
= % need to finish describing these filters
%               '': no filter; just the raw data with block-wise and
%                   trial-wise baselines subtracted
%               'Hampel': a common outlier-removing filter for pupil size
%                   data; within a sliding 2000 millisecond window,
%                   replaces the value at the center of the window with the
%                   median of the values in the window only if the original
%                   value in the center of the window is greater than 3
%                   standard deviations of the values in the window from
%                   the median of the window (feel free to play with the
%                   parameters of this filter, which are hardcoded in this
%                   function)
%               'IQR': a custom outlier-removing filter; within a sliding
%                   temporal window of 4000 milliseconds, replaces outliers
%                   with nans to exclude them from averaging into the
%                   aggregate trace, where a point at the center of the
%                   window is considered an outlier if it is +/- 3*IQR of
%                   the points in the sliding window above or below the 3rd
%                   and 1st quartiles respectively; this filter is an
%                   attempt to clean up large aberations in the pupil size
%                   data which manifest as large (and non-representative)
%                   drops in pupil size (feel free to play with the
%                   parameters of this filter, which are hardcoded in this
%                   function)
%               'SG': Savitzky-Golay smoothing filter; replaces each value
%                   with the value of the best fitting cubic polynomial
%                   within 51 frame window sliding window centered on the
%                   point to replace; this is based on the work that Scott
%                   and Maria did looking a pupil size in children with
%                   Autism
%               'windowAvg': a smoothing filter that replaces each value
%                   with the average within a sliding window of 400 ms
%                   (feel free to change the window width, which is
%                   hard-coded inside of this function); averages using the
%                   avgFunction argument, which defaults to @nanmedian
%               'BLD': a low-pass Fourier-transform-based filter; applies a
%                   trapezoidal filter to the time series data in the
%                   Fourier domain, linearly increasing attenuation
%                   between band pass and band stop frequencies of 0.5 and
%                   0.6, respectively, and differentiating the signal;
%                   thus, the resulting filtered pupil trace is is the
%                   integral of the low-pass filtered and differentiated
%                   pupil datass
%       plotWhat[optional]: choose to among the following plotting options
%           (copy the string exactly):
%               '': no plotting
%               'power': plot the low-frequency components of a Fourier
%                   power analysis (analysis from Scott)
%               'trace': plot the timecouse data of each block (as colored
%                   dotted lines) with the average across blocks
%                   superimposed as a thick blue line; verticle dotted
%                   lines mark stimulus onsets, where the white disc onsets
%                   are marked by '-.' and the black disc onset is marked
%                   by ':'; any part of the average trace which is impacted
%                   by outlier removing ('Hampel' or 'IQR') filters will be
%                   plotted in yellow
%   OUTPUT:
= % still need to decide whether to return data with baselines already subtracted
%       trialPupilSize:
%       trialBaseline:
%       blockPupilSize: 
%       blockBaseline:
%       avgPupilSizeAcrossBlocks:
%       power:


% settings
% plotWhat = 'trace'; % should be either 'power','trace', or empty
% filter = 'windowAvg'; % should be either 'Hampel','IQR','SG','windowAvg','BLD' or empty
% removeBlinks = false; % boolean
% avgFunction = @nanmedian; 

% check input
if nargin < 1
    error('AK_GABA_getPupilSizeData requires at least an asc filename (string) as an argument');
end
if nargin < 2
    directory = cd; % default to current directory
end
if nargin < 3
    removeBlinks = false; % boolean
end
if nargin < 4
    avgFunction = @nanmedian; % should be @nanmean or @nanmedian
end
if nargin < 5
    filter = ''; % should be either 'Hampel','IQR','SG','windowAvg','BLD' or empty
end
if nargin < 6
    plotWhat = ''; % should be either 'power','trace','GroupAvgTrace' or empty
end

% subject code
subj = ascFilename(1:end-6);

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
thisEdf = dir(fullfile(directory,[ascFilename(1:end-3) 'edf']));
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
        clear pSize events eventTimes blink sacc Stim1Index Stim2Index timeIndexS  % clear variables

        % extract relevant data from block struct
        time{iB} = cellfun(@str2double,block(iB).timepts(2:end,1)); % set time vector
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
            WhiteOn{iB} = eventTimes(Stim2Index)+avgError; % time of interest start times: black to white
            BlackOn{iB} = eventTimes(Stim1Index)+avgError; % white to black
        else
            WhiteOn{iB} = eventTimes(Stim2Index); % time of interest start times: black to white
            BlackOn{iB} = eventTimes(Stim1Index); % white to black
        end

        % remove far outliers from pupil data using IQR filter
        switch filter
            case 'IQR'
                % AMK IQR outliers->nan filter
                [pSize,oIdx{iB}] = AK_windowIQRfilt(pSize,4000,3); % eliminate outliers outside of the IQR by +/- 3*IQR, w/in a sliding 4000ms window  
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
                [pSize,~,~,~] = AK_BLDfilt(pSize,.5,.6,1000); % bandpass = .5; bandstop = 1; hz = 1000
        end

        % subtract out baseline
        blockBaseline(iB) = avgFunction(pSize);
        pSizeCell{iB} = pSize' - blockBaseline(iB);

        % trim pupil data timecourse to start at the first
        % BlackOn and end at last BlackOn
        pSizeCell{iB}(1:find(time{iB}==BlackOn{iB}(1))-1) = [];
        if any(strcmp(filter,{'IQR','Hampel'}))
            % trim index of points affected by filtering
            oIdx{iB}(1:find(time{iB}==BlackOn{iB}(1))-1) = [];
        end
        % trim time last so that everything is trimmed the same
        time{iB}(1:find(time{iB}==BlackOn{iB}(1))-1) = [];
        
        % index in time where events start
        timeIndexS = arrayfun(@(x) find(time{iB}==x),WhiteOn{iB});

        for TI = 1:length(WhiteOn{iB}); % cycle through times of interest: times when the white stimulus comes on  
            if removeBlinks
                % convert blink and saccade trials to nans
                if TI~=length(WhiteOn{iB}) && (~isempty(intersect(blink,(WhiteOn{iB}(TI)-500:WhiteOn{iB}(TI+1)))) || ~isempty(intersect(sacc,(WhiteOn{iB}(TI)-500:WhiteOn{iB}(TI+1))))) 
                    pSizeCell{iB}(timeIndexS(TI):timeIndexS(TI+1)-1) = nan;
                elseif ~isempty(intersect(blink,(WhiteOn{iB}(TI)-500:time{iB}(end)))) || ~isempty(intersect(sacc,(WhiteOn{iB}(TI)-500:time{iB}(end))))
                    pSizeCell{iB}(timeIndexS(TI):end) = nan;
                end
            end
            % relative pupil size trace for each trial
            try
                trialBaseline(TI+10*(iB-1)) = pSizeCell{iB}(timeIndexS(TI)-250);
                trialPupilSize(TI+10*(iB-1),:) = pSizeCell{iB}(timeIndexS(TI):timeIndexS(TI)+1999) - trialBaseline(TI+10*(iB-1));
            catch
                if ~isempty(pSizeCell{iB}(timeIndexS(TI):end))
                    trialBaseline(TI+10*(iB-1)) = pSizeCell{iB}(timeIndexS(TI)-250);
                    trialPupilSize(TI+10*(iB-1),:) = [pSizeCell{iB}(timeIndexS(TI):end) nan(1,2000-length(pSizeCell{iB}(timeIndexS(TI):end)))] - trialBaseline(TI+10*(iB-1));
                else
                    trialBaseline(TI+10*(iB-1)) = nan;
                    trialPupilSize(TI+10*(iB-1),:) = nan(1,2000);
                end
            end
        end
    end
end

% set up figures for traces per subject
if strcmp(plotWhat,'trace')
    figure; hold on;
    blockColors = {'r','g','m','c'};
end

% nan pad to make each block the same length for concatenation
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
        plot(1:length(pSizeCell{iB}),pSizeCell{iB},['--' blockColors{iB}]); % no extra filtering
        % plot event code loci
        for iW = 1:length(WhiteOn{iB})
            plot([WhiteOn{iB}(iW)-time{iB}(1) WhiteOn{iB}(iW)-time{iB}(1)],[min(cell2mat(cellfun(@min,pSizeCell,'UniformOutput',false))) max(cell2mat(cellfun(@max,pSizeCell,'UniformOutput',false)))],'-.k');
        end
        for iBl = 1:length(BlackOn{iB})
            plot([BlackOn{iB}(iBl)-time{iB}(1) BlackOn{iB}(iBl)-time{iB}(1)],[min(cell2mat(cellfun(@min,pSizeCell,'UniformOutput',false))) max(cell2mat(cellfun(@max,pSizeCell,'UniformOutput',false)))],':k');
        end
    end
end

% cat
clear pupilData 
blockPupilSize = vertcat(pSizeCell{:});

% avg across blocks
avgPupilSizeAcrossBlocks = avgFunction(blockPupilSize);

if strcmp(plotWhat,'trace')
    % create vector of timepoints since stim onset for x-axis
    plotTime = 1:length(avgPupilSizeAcrossBlocks);                              

    % plot avg size at time point
    plot(plotTime,avgPupilSizeAcrossBlocks,'-b','LineWidth',2);
    if any(strcmp(filter,{'IQR','Hampel'}))
        % show where the outlier detection filter has affected the average
        clear hampelIdx hampelTrace
        hampelIdx = sum(vertcat(oIdx{:})) > 0;
        hampelTrace = avgPupilSizeAcrossBlocks;
        hampelTrace(~hampelIdx) = nan;
        plot(plotTime,hampelTrace,'-y','LineWidth',2);
    end
    aH = gca; 

    % manipulate axis properties
    title = [subj ' Avg Pupil Size Across Blocks'];
    if ~isempty(filter)
        title = [title ': ' filter ' filter'];
    end
    set(aH,'Title',text('String',title));
    xlabel('Time (ms)');
    ylabel('Average Pupil Size');
    axis([min(plotTime) max(plotTime) min(avgPupilSizeAcrossBlocks)-200 max(avgPupilSizeAcrossBlocks)+200]);
    set(aH,'XTick',1:2000:max(plotTime));
end

% Scott's Fourier analysis
clear y f t p power
[y,f,t,p] = spectrogram(nanmean(blockPupilSize),1000,0,128,1000);
power = 10*log10(abs(p));

if strcmp(plotWhat,'power')
    figure
    surf(t,f,10*log10(abs(p)),'EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90);
    xlabel('Time');
    ylabel('Frequency (Hz)');
    axis([.5 20 0 50])
    title = [subj ' Power Analysis Across Blocks'];
    set(gca,'Title',text('String',title));

    figure
    hold on
    plot(power(1,:))
    xlabel('Time');
    ylabel('Power @ low Hz');
    title = [subj ' Power At Low Frequencies'];
    set(gca,'Title',text('String',title));
end

% invert baseline arrays so that dimensions match pupilSize arrays
trialBaseline = trialBaseline';
blockBaseline = blockBaseline';

% save data matrix for Scott
% if ~isempty(filter)
%     save(fullfile(directory,[subj '_PupilSize_' filter 'filter.mat']),'pupilData');
% else
%     save(fullfile(directory,[subj '_PupilSize.mat']),'pupilData');
% end

end

