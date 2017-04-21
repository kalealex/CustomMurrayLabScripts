% Psychophysics_Pupil_Fatigue_Lorazepam
% Murray Lab 2016
% created by AMK 12/20/16

%% settings

loadMetrics = false; % load metrics (1) or recalculate and save (0)
 removeBlinks = true; % interpolate across blinks? yes (1) or no (0)
 chunkSize = 82000; % milliseconds
 experiment = 'o'; % motion ('m'), contrast ('o'), orientation ('o')
classifyDifference = true; % classify differences (1) between session stats or classify sessions individually (0)
weightMeans = true; % weight means by block length? yes(1) or no (0)
blindLast5Subj = false; % training set for LDA and cross-validation accuracies based on first 10 subjects (true) or all subjects (false, no test set)
 if blindLast5Subj
     trainingN = 10;
 else 
     trainingN = 15;
 end

% this will use all possible combinations of metrics from all experiment types to make a table of crossvalidation accuracies and predictions
makePredictionTable = true;  
 if makePredictionTable
     loadMetrics = true;
     experiment = 'NA';
     classifyDifference = true;
 end

%% establish directory by subjects

top_dir = 'L:\MurrayLab\Lorazepam\Data'; % set up base directory
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab\MotionStairPsychophysics_EyetrackingResponse\Lorazepam'; % set directory for saving data and figures
subj_dirs = {'L001','L002','L003','L004','L005','L006','L007','L008','L009','L010','L011','L012','L013','L014','L015'};
% subj_dirs = {'L001','L002','L003','L004','L005','L006','L007','L008','L009','L010','L011','L012','L013'};
session_dirs = {'Session 1','Session 2'};

%% parse .asc files and calculate metrics on fatigue from pupil size data

if ~loadMetrics
    % set expected number of events
    if strcmp(experiment,'m')
        eventsExpected = 572;
    elseif strcmp(experiment,'c')
        eventsExpected = 482;
    elseif strcmp(experiment,'o')
        eventsExpected = 2;
    else
        error('Invalid experiment (see settings)')
    end
    for iSubj = 1:length(subj_dirs); % cycle through subjects
        subj_dir = fullfile(top_dir,subj_dirs{iSubj});
        if exist(subj_dir,'dir')==0
            disp('could not find subject folder') % message
        elseif exist(subj_dir,'dir')==7
            % create waitbar
            wH = waitbar(0,['Creating Directory for ' subj_dirs{iSubj}]); 
            for iSess = 1:length(session_dirs)
                % generate directory in use
                use_dir = fullfile(subj_dir,session_dirs{iSess}); 
                asc_list = dir(fullfile(use_dir,['*' experiment '*.asc'])); % generate structure containing information about pupil asc files
                if ~isempty(asc_list)
                    for iA = 1:length(asc_list); % cycle through asc files
                        % waitbar
                        wait = (iSess+length(session_dirs)*(iSubj-1))/(length(subj_dirs)*length(session_dirs)+1); % set waitbar progress
                        waitbar(wait,wH,['processing ' asc_list(iA).name(1:end-4)]); % update waitbar
                        % load data
                        clear block mat_name thisBlock
                        mat_name = [fullfile(use_dir,asc_list(iA).name(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
                        if ~exist(mat_name,'file') % check whether or not function output has already been saved
                            disp(['parsing ' fullfile(use_dir,asc_list(iA).name)]) % message
                            block = AK_GABAinASD_Psychophysics_ascfileParse(fullfile(use_dir,asc_list(iA).name)); % parse asc file
                            save(mat_name,'block'); % save function output
                        else
                            disp(['loading ' mat_name]) % message
                            load(mat_name,'block') % load file containing timepts array
                        end
                        % store metadata
                        ascFilename{iSubj,iSess} = asc_list(iA).name; % document name of file of interest for saving figures
                        nBlocks(iSubj,iSess) = length(block); % document number of blocks of trials in file
                        % set up block counter
                        thisBlock = 0;
                        for iB = 1:length(block); % cycle through blocks of trials
                            % only proceed if there are a full set events
                            % and if time is continuous
                            %and if duration is longer than 2 * chunkSize 
                            if length(block(iB).events(:,1)) == eventsExpected &&...
                                    ~any((cellfun(@str2double,block(iB).timepts(3:end,1))-cellfun(@str2double,block(iB).timepts(2:end-1,1)))~=1) &&...
                                    length(block(iB).timepts(2:end,1)) >= 2 * chunkSize
                                clear time pSize pSizeFilt events eventTimes blink Stim1Index Stim2Index timeIndexS  % clear variables
                                
                                % increase block count
                                thisBlock = thisBlock + 1;
                                
                                % extract relevant data from block struct
                                time = cellfun(@str2double,block(iB).timepts(2:end,1)); % set time vector
                                pSize = cellfun(@str2double,block(iB).timepts(2:end,4)); % set pupil size vector

                                events = block(iB).events(2:end,2); % store list of events as cell array of strings
                                eventTimes = cellfun(@str2double,block(iB).events(2:end,1)); % set vector of event times
                                events(cellfun(@isempty,events)) = {'no event'}; % fill in empty events
                                
                                % store block length for weighting avgs
                                if isempty(length(time))
                                    blockLength(iSubj,iSess,thisBlock) = 0;
                                else
                                    blockLength(iSubj,iSess,thisBlock) = length(time);
                                end
                                
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
                                    % remove blinks and saccades from signal by interpolating
                                    pSize(ismember(time,blink)) = nan; pSize(ismember(time,sacc)) = nan;
                                    pSize = AK_interpolateAcrossNans(pSize);
                                end

                                % calculate miosis and hippus metrics
                                [PUI(iSubj,iSess,thisBlock), CM(iSubj,iSess,thisBlock), PDRmin(iSubj,iSess,thisBlock)] = AK_PupilSize_FatigueMetrics_Miosis(pSize,chunkSize); % miosis metrics for 82 second chunks
                                PowerLowFreq(iSubj,iSess,thisBlock) = AK_PupilSize_FatigueMetrics_Hippus(pSize,1000); % hippus metric for sample rate of 1000 hz


    %                             % remove far outliers from pupil data using IQR filter
    %                             switch filter
    %                                 case 'IQR'
    %                                     % AMK IQR outliers->nan filter
    %                                     [pSize,oIdx{thisBlock}] = AK_windowIQRfilt(pSize,4000,3); % eliminate outliers outside of the IQR by +/- 1.5*IQR, w/in a sliding 2000ms window  
    %                                 case 'Hampel'
    %                                     % Hampel outliers->windowMedian filter
    %                                     [pSize,oIdx{thisBlock}] = hampel(pSize,2000,3);
    %                                 case 'SG'
    %                                     % apply Savitzky-Golay smoothing filter; local cubic modeling w/in 51 frame window 
    %                                     pSize = sgolayfilt(pSize,3,51); 
    %                                       % modified frame rate from Maria's code, where param = 13:
    %                                       % 51=13*(1000hz/250hz)-1 (frame rate must be odd)
    %                                 case 'windowAvg'
    %                                     % AMK window average smoothing filter
    %                                     pSize = AK_windowAvgFilt(pSize,400,avgFunction);
    %                                 case 'BLD'
    %                                     % AMK band-limited differentiator filter (Low-Pass)
    %                                     [~,pSize,~] = AK_BLDfilt(pSize,.5,.6,1000); % bandpass = .5; bandstop = 1; hz = 1000
    %                             end


                            end
                        end
                    end
                    wait = 1; % set waitbar progress
                    waitbar(wait,wH,'Finishing'); % update waitbar
                else
                    disp('could not find .asc files') % message
                end
            end
            close(wH);
        end
    end
else
    if strcmp(experiment,'m')
        load(fullfile(home_dir,'ssMotionLorazepam_PupilFatigueMetrics.mat'));
    elseif strcmp(experiment,'c')
        load(fullfile(home_dir,'ContrastLorazepam_PupilFatigueMetrics.mat'));
    elseif strcmp(experiment,'o')
        load(fullfile(home_dir,'OrientationLorazepam_PupilFatigueMetrics.mat'));
    end
end

%% guess which session was drug vs placebo based on pupil size

if ~makePredictionTable
    % avg across blocks
    if ~weightMeans
        mPUI = mean(PUI,3);
        mCM = mean(CM,3);
        mPDRmin = mean(PDRmin,3);
        mPowerLowFreq = mean(PowerLowFreq,3);
    else
        mPUI = wmean(PUI,blockLength,3);
        mCM = wmean(CM,blockLength,3);
        mPDRmin = wmean(PDRmin,blockLength,3);
        mPowerLowFreq = wmean(PowerLowFreq,blockLength,3);
    end

    % make a prediction of which session was which for each participant (guess
    % == 1 means first session was drug; guess == 0 means second session was
    % drug)
    for iS = 1:length(subj_dirs)
       PUIguess(iS) = mPUI(iS,1) > mPUI(iS,2); % session with drug should have greater pupillary unrest index
       CMguess(iS) = mCM(iS,1) > mCM(iS,2); % session with drug should have greater cumulative miosis
       PDRminGuess(iS) = mPDRmin(iS,1) < mPDRmin(iS,2); % session with drug should have smaller minimum pupillary dilation ratio 
       PowerLowFreqGuess(iS) = mPowerLowFreq(iS,1) > mPowerLowFreq(iS,2); % session with drug should have greater power at low frequency oscillations
    end

    guessKey = 'The variables with guess in their name are logical arrays of predictions per subject about whether drug (1) or placebo (0) was administered on Session 1';
end

%% save metrics

if ~loadMetrics
    if strcmp(experiment,'m')
        save(fullfile(home_dir,'ssMotionLorazepam_PupilFatigueMetrics.mat'),'PUI','CM','PDRmin','PowerLowFreq','blockLength','PUIguess','CMguess','PDRminGuess','PowerLowFreqGuess','guessKey');
    elseif strcmp(experiment,'c')
        save(fullfile(home_dir,'ContrastLorazepam_PupilFatigueMetrics.mat'),'PUI','CM','PDRmin','PowerLowFreq','blockLength','PUIguess','CMguess','PDRminGuess','PowerLowFreqGuess','guessKey');
    elseif strcmp(experiment,'o')
        save(fullfile(home_dir,'OrientationLorazepam_PupilFatigueMetrics.mat'),'PUI','CM','PDRmin','PowerLowFreq','blockLength','PUIguess','CMguess','PDRminGuess','PowerLowFreqGuess','guessKey');
    end
end

%% load actual data to compare guess with correct answer

% % load true session info
% load('L:\MurrayLab\Lorazepam\analysis_code\blind_order.mat');
% 
% % index of subjective inebriation 
% hardToTellIdx = [1 4 10];
% inebriationIdx = [3 5 9];
% 
% % organize statistics by drug vs placebo
% for iS = 1:size(blind_order,2)
%     drug(iS).PUI = mPUI(iS,blind_order(1,iS));
%     drug(iS).CM = mCM(iS,blind_order(1,iS));
%     drug(iS).PDRmin = mPDRmin(iS,blind_order(1,iS));
%     drug(iS).PowerLowFreq = mPowerLowFreq(iS,blind_order(1,iS));
%     placebo(iS).PUI = mPUI(iS,blind_order(2,iS));
%     placebo(iS).CM = mCM(iS,blind_order(2,iS));
%     placebo(iS).PDRmin = mPDRmin(iS,blind_order(2,iS));
%     placebo(iS).PowerLowFreq = mPowerLowFreq(iS,blind_order(2,iS));
% end
% 
% % create matrices for each summary statistic (this is how Scott and I looked at the data in our meeting)
% matPUI = [drug.PUI; placebo.PUI]';
% matCM = [drug.CM; placebo.CM]';
% matPDRmin = [drug.PDRmin; placebo.PDRmin]';
% matPowerLowFreq = [drug.PowerLowFreq; placebo.PowerLowFreq]';
% 
% % differences
% diffPUI = (matPUI(:,1) - matPUI(:,2))./matPUI(:,2);
% diffCM = (matCM(:,1) - matCM(:,2))./matCM(:,2);
% diffPDRmin = (matPDRmin(:,1) - matPDRmin(:,2))./matPDRmin(:,2);
% diffPowerLowFreq = (matPowerLowFreq(:,1) - matPowerLowFreq(:,2))./matPowerLowFreq(:,2);

%% use linear discriminant analysis to predict the drug vs placebo sessions

% load true session info
if ~exist('blind_order','var')
    load('L:\MurrayLab\Lorazepam\analysis_code\blind_order.mat');
end

if ~makePredictionTable
    clear trainingSet groupLabel testSet

    % organize data for input to classify function:
    if ~classifyDifference
        for iS = 1:size(mPUI,1)
            if iS <= trainingN
                % use first ten participants' data as training set
                trainingSet(2*(iS-1)+1:2*(iS-1)+2,1) = mPUI(iS,:)'; % PUI
                trainingSet(2*(iS-1)+1:2*(iS-1)+2,2) = mCM(iS,:)'; % CM
                trainingSet(2*(iS-1)+1:2*(iS-1)+2,3) = mPDRmin(iS,:)'; % PDRmin
                trainingSet(2*(iS-1)+1:2*(iS-1)+2,4) = mPowerLowFreq(iS,:)'; % Power
                % label each point with a drug or placebo
                if blind_order(1,iS)==1 % first session for this subbject was drug
                    groupLabel(2*(iS-1)+1:2*(iS-1)+2,1) = {'drug';'placebo'};
                elseif blind_order(1,iS)==2 % second session for this subbject was drug
                    groupLabel(2*(iS-1)+1:2*(iS-1)+2,1) = {'placebo';'drug'};
                end
            else
                % use other participants' as a test set
                testSet(2*(iS-11)+1:2*(iS-11)+2,1) = mPUI(iS,:)'; % PUI
                testSet(2*(iS-11)+1:2*(iS-11)+2,2) = mCM(iS,:)'; % CM
                testSet(2*(iS-11)+1:2*(iS-11)+2,3) = mPDRmin(iS,:)'; % PDRmin
                testSet(2*(iS-11)+1:2*(iS-11)+2,4) = mPowerLowFreq(iS,:)'; % Power
            end
        end
    else
        for iS = 1:size(mPUI,1)
            if iS <= 10
                % use first ten participants' data as training set
                trainingSet(iS,1) = mPUI(iS,1)-mPUI(iS,2); % PUI
                trainingSet(iS,2) = mCM(iS,1)-mCM(iS,2); % CM
                trainingSet(iS,3) = mPDRmin(iS,1)-mPDRmin(iS,2); % PDRmin
                trainingSet(iS,4) = mPowerLowFreq(iS,1)-mPowerLowFreq(iS,2); % Power
                % label each point with a drug or placebo
                if blind_order(1,iS)==1 % first session for this subbject was drug
                    groupLabel(iS,1) = {'drug'};
                elseif blind_order(1,iS)==2 % first session for this subbject was placebo
                    groupLabel(iS,1) = {'placebo'};
                end
            else
                % use other participants' as a test set
                testSet(iS-trainingN,1) = mPUI(iS,1)-mPUI(iS,2); % PUI
                testSet(iS-trainingN,2) = mCM(iS,1)-mCM(iS,2); % CM
                testSet(iS-trainingN,3) = mPDRmin(iS,1)-mPDRmin(iS,2); % PDRmin
                testSet(iS-trainingN,4) = mPowerLowFreq(iS,1)-mPowerLowFreq(iS,2); % Power
            end
        end
    end
    
    % leave-one-out cross validation
    accuracy = AK_leaveOneOutCrossValidation(trainingSet,groupLabel);
    % run linear discriminant analysis
    [class,~,~,~,coeff] = classify(testSet,trainingSet,groupLabel,'linear');
elseif makePredictionTable    
    % designate experiments to use
    useExperiments = {'m','c','o','all'};
    % cycle through experiments loading data and creating a set of predictions for each
    for iE = 1:length(useExperiments)
        % load metrics and get avgs
        clear PUI CM PDRmin PowerLowFreq mPUI mCM mPDRmin mPowerLowFreq tableRow accuracy class vars tempTable
        if strcmp(useExperiments(iE),'m')
            load(fullfile(home_dir,'ssMotionLorazepam_PupilFatigueMetrics.mat'));
            % avg across blocks
            if ~weightMeans
                mPUI = mean(PUI,3);
                mCM = mean(CM,3);
                mPDRmin = mean(PDRmin,3);
                mPowerLowFreq = mean(PowerLowFreq,3);
            else
                mPUI = wmean(PUI,blockLength,3);
                mCM = wmean(CM,blockLength,3);
                mPDRmin = wmean(PDRmin,blockLength,3);
                mPowerLowFreq = wmean(PowerLowFreq,blockLength,3);
            end
        elseif strcmp(useExperiments(iE),'c')
            load(fullfile(home_dir,'ContrastLorazepam_PupilFatigueMetrics.mat'));
            % avg across blocks
            if ~weightMeans
                mPUI = mean(PUI,3);
                mCM = mean(CM,3);
                mPDRmin = mean(PDRmin,3);
                mPowerLowFreq = mean(PowerLowFreq,3);
            else
                mPUI = wmean(PUI,blockLength,3);
                mCM = wmean(CM,blockLength,3);
                mPDRmin = wmean(PDRmin,blockLength,3);
                mPowerLowFreq = wmean(PowerLowFreq,blockLength,3);
            end
        elseif strcmp(useExperiments(iE),'o')
            load(fullfile(home_dir,'OrientationLorazepam_PupilFatigueMetrics.mat'));
            % avg across blocks
            if ~weightMeans
                mPUI = mean(PUI,3);
                mCM = mean(CM,3);
                mPDRmin = mean(PDRmin,3);
                mPowerLowFreq = mean(PowerLowFreq,3);
            else
                mPUI = wmean(PUI,blockLength,3);
                mCM = wmean(CM,blockLength,3);
                mPDRmin = wmean(PDRmin,blockLength,3);
                mPowerLowFreq = wmean(PowerLowFreq,blockLength,3);
            end
        elseif strcmp(useExperiments(iE),'all')
            % load all and average all together:
            load(fullfile(home_dir,'ssMotionLorazepam_PupilFatigueMetrics.mat')); % motion
            % avg across blocks
            if ~weightMeans
                motionPUI = mean(PUI,3);
                motionCM = mean(CM,3);
                motionPDRmin = mean(PDRmin,3);
                motionPowerLowFreq = mean(PowerLowFreq,3);
            else
                motionPUI = wmean(PUI,blockLength,3);
                motionCM = wmean(CM,blockLength,3);
                motionPDRmin = wmean(PDRmin,blockLength,3);
                motionPowerLowFreq = wmean(PowerLowFreq,blockLength,3);
                motionBlockLengths = mean(blockLength,3);
            end
            clear PUI CM PDRmin PowerLowFreq
            load(fullfile(home_dir,'ContrastLorazepam_PupilFatigueMetrics.mat')); % contrast
            % avg across blocks
            if ~weightMeans
                contrastPUI = mean(PUI,3);
                contrastCM = mean(CM,3);
                contrastPDRmin = mean(PDRmin,3);
                contrastPowerLowFreq = mean(PowerLowFreq,3);
            else
                contrastPUI = wmean(PUI,blockLength,3);
                contrastCM = wmean(CM,blockLength,3);
                contrastPDRmin = wmean(PDRmin,blockLength,3);
                contrastPowerLowFreq = wmean(PowerLowFreq,blockLength,3);
                contrastBlockLengths = mean(blockLength,3);
            end
            clear PUI CM PDRmin PowerLowFreq
            load(fullfile(home_dir,'OrientationLorazepam_PupilFatigueMetrics.mat'));
            % avg across blocks
            if ~weightMeans
                orientPUI = mean(PUI,3);
                orientCM = mean(CM,3);
                orientPDRmin = mean(PDRmin,3);
                orientPowerLowFreq = mean(PowerLowFreq,3);
            else
                orientPUI = wmean(PUI,blockLength,3);
                orientCM = wmean(CM,blockLength,3);
                orientPDRmin = wmean(PDRmin,blockLength,3);
                orientPowerLowFreq = wmean(PowerLowFreq,blockLength,3);
                orientBlockLengths = mean(blockLength,3);
            end
            % grand avgs
            if ~weightMeans
                mPUI = mean(cat(3,motionPUI,contrastPUI,orientPUI),3);
                mCM = mean(cat(3,motionCM,contrastCM,orientCM),3);
                mPDRmin = mean(cat(3,motionPDRmin,contrastPDRmin,orientPDRmin),3);
                mPowerLowFreq = mean(cat(3,motionPowerLowFreq,contrastPowerLowFreq,orientPowerLowFreq),3);
            else
                mPUI = wmean(cat(3,motionPUI,contrastPUI,orientPUI),cat(3,motionBlockLengths,contrastBlockLengths,orientBlockLengths),3);
                mCM = wmean(cat(3,motionCM,contrastCM,orientCM),cat(3,motionBlockLengths,contrastBlockLengths,orientBlockLengths),3);
                mPDRmin = wmean(cat(3,motionPDRmin,contrastPDRmin,orientPDRmin),cat(3,motionBlockLengths,contrastBlockLengths,orientBlockLengths),3);
                mPowerLowFreq = wmean(cat(3,motionPowerLowFreq,contrastPowerLowFreq,orientPowerLowFreq),cat(3,motionBlockLengths,contrastBlockLengths,orientBlockLengths),3);
            end
        end

        clear trainingSet groupLabel testSet

        % organize data for input to classify function:
        if ~classifyDifference
            for iS = 1:size(mPUI,1)
                if iS <= trainingN
                    % use first ten participants' data as training set
                    trainingSet(2*(iS-1)+1:2*(iS-1)+2,1) = mPUI(iS,:)'; % PUI
                    trainingSet(2*(iS-1)+1:2*(iS-1)+2,2) = mCM(iS,:)'; % CM
                    trainingSet(2*(iS-1)+1:2*(iS-1)+2,3) = mPDRmin(iS,:)'; % PDRmin
                    trainingSet(2*(iS-1)+1:2*(iS-1)+2,4) = mPowerLowFreq(iS,:)'; % Power
                    % label each point with a drug or placebo
                    if blind_order(1,iS)==1 % first session for this subbject was drug
                        groupLabel(2*(iS-1)+1:2*(iS-1)+2,1) = {'drug';'placebo'};
                    elseif blind_order(1,iS)==2 % second session for this subbject was drug
                        groupLabel(2*(iS-1)+1:2*(iS-1)+2,1) = {'placebo';'drug'};
                    end
                else
                    % use other participants' as a test set
                    testSet(2*(iS-11)+1:2*(iS-11)+2,1) = mPUI(iS,:)'; % PUI
                    testSet(2*(iS-11)+1:2*(iS-11)+2,2) = mCM(iS,:)'; % CM
                    testSet(2*(iS-11)+1:2*(iS-11)+2,3) = mPDRmin(iS,:)'; % PDRmin
                    testSet(2*(iS-11)+1:2*(iS-11)+2,4) = mPowerLowFreq(iS,:)'; % Power
                end
            end
        else
            for iS = 1:size(mPUI,1)
                if iS <= trainingN
                    % use first ten participants' data as training set
                    trainingSet(iS,1) = mPUI(iS,1)-mPUI(iS,2); % PUI
                    trainingSet(iS,2) = mCM(iS,1)-mCM(iS,2); % CM
                    trainingSet(iS,3) = mPDRmin(iS,1)-mPDRmin(iS,2); % PDRmin
                    trainingSet(iS,4) = mPowerLowFreq(iS,1)-mPowerLowFreq(iS,2); % Power
                    % label each point with a drug or placebo
                    if blind_order(1,iS)==1 % first session for this subbject was drug
                        groupLabel(iS,1) = {'drug'};
                    elseif blind_order(1,iS)==2 % first session for this subbject was placebo
                        groupLabel(iS,1) = {'placebo'};
                    end
                else
                    % use other participants' as a test set
                    testSet(iS-trainingN,1) = mPUI(iS,1)-mPUI(iS,2); % PUI
                    testSet(iS-trainingN,2) = mCM(iS,1)-mCM(iS,2); % CM
                    testSet(iS-trainingN,3) = mPDRmin(iS,1)-mPDRmin(iS,2); % PDRmin
                    testSet(iS-trainingN,4) = mPowerLowFreq(iS,1)-mPowerLowFreq(iS,2); % Power
                end
            end
        end

        % create table of all possible combinations of metrics and corresponding accuracies and predicitions:
        metricLabels = {'PUI','CM','PDRmin','Power'}; % same order as above
        % all possible combinations of indices
        choose4 = nchoosek([1 2 3 4],4); 
        choose3 = nchoosek([1 2 3 4],3);
        choose2 = nchoosek([1 2 3 4],2);
        choose1 = nchoosek([1 2 3 4],1);
        % set up counter
        iRow = 1;
        % run all combinations of 4
        for iC4 = 1:size(choose4,1)
            % cat row labels
            tableRows{iRow} = strjoin(metricLabels(choose4(iC4,:)),', ');
            % leave-one-out cross validation
            accuracy(iRow) = AK_leaveOneOutCrossValidation(trainingSet(:,choose4(iC4,:)),groupLabel);
            % run linear discriminant analysis
            try
                [class{iRow},~,~,~,~] = classify(testSet(:,choose4(iC4,:)),trainingSet(:,choose4(iC4,:)),groupLabel,'linear');
            catch
                class{iRow} = 'LDA failed';
            end
            % increment counter
            iRow = iRow + 1;
        end
        % run all combinations of 3
        for iC3 = 1:size(choose3,1)
            % cat row labels
            tableRows{iRow} = strjoin(metricLabels(choose3(iC3,:)),', ');
            % leave-one-out cross validation
            accuracy(iRow) = AK_leaveOneOutCrossValidation(trainingSet(:,choose3(iC3,:)),groupLabel);
            % run linear discriminant analysis
            try
                [class{iRow},~,~,~,~] = classify(testSet(:,choose3(iC3,:)),trainingSet(:,choose3(iC3,:)),groupLabel,'linear');
            catch
                class{iRow} = 'LDA failed';
            end
            % increment counter
            iRow = iRow + 1;
        end
        % run all combinations of 2
        for iC2 = 1:size(choose2,1)
            % cat row labels
            tableRows{iRow} = strjoin(metricLabels(choose2(iC2,:)),', ');
            % leave-one-out cross validation
            accuracy(iRow) = AK_leaveOneOutCrossValidation(trainingSet(:,choose2(iC2,:)),groupLabel);
            % run linear discriminant analysis
            try
                [class{iRow},~,~,~,~] = classify(testSet(:,choose2(iC2,:)),trainingSet(:,choose2(iC2,:)),groupLabel,'linear');
            catch
                class{iRow} = 'LDA failed';
            end
            % increment counter
            iRow = iRow + 1;
        end
        % run all combinations of 1
        for iC1 = 1:size(choose1,1)
            % cat row labels
            tableRows{iRow} = strjoin(metricLabels(choose1(iC1,:)),', ');
            % leave-one-out cross validation
            accuracy(iRow) = AK_leaveOneOutCrossValidation(trainingSet(:,choose1(iC1,:)),groupLabel);
            % run linear discriminant analysis
            try
                [class{iRow},~,~,~,~] = classify(testSet(:,choose1(iC1,:)),trainingSet(:,choose1(iC1,:)),groupLabel,'linear');
            catch
                class{iRow} = 'LDA failed';
            end
            % increment counter
            iRow = iRow + 1;
        end
        % make table/ append to table
        vars = {[useExperiments{iE} 'Accuracy'],[useExperiments{iE} 'Classifications']}; % column labels
        if iE == 1 % create
            lzPupilMetricsTable = table(accuracy',class','RowNames',tableRows','VariableNames',vars);
        else % append
            tempTable = table(accuracy',class','RowNames',tableRows','VariableNames',vars);
            lzPupilMetricsTable = [lzPupilMetricsTable tempTable];
        end
    end
    % save table
    if blindLast5Subj && weightMeans
        save(fullfile(home_dir,'Lorazepam_weightedPupilFatigueMetrics_blindedPredictionTable.mat'),'lzPupilMetricsTable');
    elseif blindLast5Subj && ~weightMeans
        save(fullfile(home_dir,'Lorazepam_PupilFatigueMetrics_blindedPredictionTable.mat'),'lzPupilMetricsTable');
    elseif ~blindLast5Subj && weightMeans
        save(fullfile(home_dir,'Lorazepam_weightedPupilFatigueMetrics_unblindedPredictionTable.mat'),'lzPupilMetricsTable');
    else
        save(fullfile(home_dir,'Lorazepam_PupilFatigueMetrics_unblindedPredictionTable.mat'),'lzPupilMetricsTable');
    end
end

