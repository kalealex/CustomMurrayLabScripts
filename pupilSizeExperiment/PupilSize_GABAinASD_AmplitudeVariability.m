% PupilSize_GABAinASD_AmplitudeVariability
% MurrayLab_2017
% created by AMK 20170511

%% designate root directory

root_dir = 'L:\MurrayLab\ASD\Data';

%% load previously saved trial exclusions?

usePreviouslySavedExclusion = true;
previousSaveFile = 'pupilAmpData2.mat';

if usePreviouslySavedExclusion
    load(fullfile('C:\Users\Alex Kale\Documents\MATLAB\CustomMurrayLabScripts\pupilSizeExperiment\amplitudeAnalysisData',previousSaveFile));
    oldSubjects = subjects;
    oldExcludeTrial = excludeTrial;
end

%% find subjects to include in analysis

% for GABA in ASD
subjects = dir(fullfile(root_dir,'G*')); % structure of all matching items in directory
subjects = {subjects(:).name}; % cell array of subject codes/ folder names
subjects = [subjects {'KG122'}]; % add kiddo
subjects(cellfun(@(x) exist(x,'dir')==0,cellfun(@(x) fullfile(root_dir,x),subjects,'UniformOutput',false))) = []; % remove elements which are not valid directories

if usePreviouslySavedExclusion
    % use old exclusion index but add in any new subjects
    newSubjects = setdiff(subjects,oldSubjects);
end

%% settings for trial exclusion and amplitude calculation

velocityCriterion = 700; % percent pupil size change within 1 ms 
hz = 1000; % sampling rate of eye tracker
dt = 1/hz; % time differential for calculating velocity trace

ampWindowHalfWidth = 25; % half of window width for averaging to find amplitude

%% get amplitude variability for each subject

% preallocate
excludeTrial = false(length(subjects),40); % exclusion log per trial per subject
amp = nan(length(subjects),40); % amp per trial per subject
timeMax = nan(length(subjects),40); % time of peak pupil size per trial per subject
timeMin = nan(length(subjects),40); % time of trough pupil size per trial per subject
ampSD = nan(length(subjects),1); % standard deviations of amp per subject
timeMaxSD = nan(length(subjects),1); % standard deviations of time of peak pupil size per subject
timeMinSD = nan(length(subjects),1); % standard deviations of time of trough pupil size per subject
% trimmedAmp = nan(length(subjects),40); % amp per trial per subject, trimmed 20% of the highest amplitude
% trimmedAmpSD = nan(length(subjects),1); % standard deviations of amp per subject, trimmed 20% of the highest amplitude
% allTrials = [];

% empty figure so that old trials are not confused for the next trial to
% judge
figure;

for iSubj = 1:length(subjects)
    % navigate to sub-directory and check that it exists
    data_dir = fullfile(root_dir,subjects{iSubj},'Psychophysics');
    if exist(data_dir,'dir') ~= 0
        cd(data_dir);
        % get .asc files that correspond to pupil data
        asc_list = dir([subjects{iSubj} 'p*.asc']);
        if ~isempty(asc_list)
            asc_list = AK_sortStruct(asc_list,1,[subjects{iSubj} 'p']);
            disp(asc_list.name); % test
            for iAsc = 1:length(asc_list)
                if exist(asc_list(iAsc).name,'file')
                    % parse data
                    [~, data] = AK_GABA_getPupilSizeData(asc_list(iAsc).name,cd,true,@nanmedian,'Hampel',''); % get data, remove blinks, use median averaging and Hampel filter
                    % document amplitudes for this subject:
                    % Scott idea for amplitude (peak median - trough
                    % median), where peak and trough are defined
                    % trial-by-trial; exclude trials using velocity
                    % threshold and inference by eye
                    for iTrial = 1:length(data.pupilTracePerTrial(:,1))
                        clear maxIdx minIdx maxWindow minWindow velocity                        
                        % get temporal windows for averaging to calculate amplitude
                        [~, maxIdx] = max(data.pupilTracePerTrial(iTrial,1:1000)); % find max/min value w/in first second of trial
                        [~, minIdx] = min(data.pupilTracePerTrial(iTrial,1:1000));
                        if maxIdx - ampWindowHalfWidth < 1
                            maxWindow = 1:(maxIdx + 2*ampWindowHalfWidth - (length(1:maxIdx)-1));
                        else
                            maxWindow = (maxIdx-ampWindowHalfWidth):(maxIdx+ampWindowHalfWidth);
                        end
                        if minIdx - ampWindowHalfWidth < 1
                            minWindow = 1:(minIdx + 2*ampWindowHalfWidth - (length(1:minIdx)-1));
                        else
                            minWindow = (minIdx-ampWindowHalfWidth):(minIdx+ampWindowHalfWidth);
                        end
                        % check whether this trial is worth including
                        if usePreviouslySavedExclusion && ~any(strcmp(newSubjects,subjects{iSubj}))
                            % use previously saved exclusions only for
                            % subjects whose data has been examined
                            % previously 
                            excludeTrial(iSubj,iTrial) = oldExcludeTrial(strcmp(subjects{iSubj},oldSubjects),iTrial);
                        else
                            % examine trace if velocity threshold is
                            % exceeded
                            velocity = diff(data.pupilTracePerTrial(iTrial,:),1)/dt;
    %                         [~,velocity,~,~] = AK_BLDfilt(data.pupilTracePerTrial(iTrial,:)',.5,.6,hz);
                            if any(abs(velocity) > velocityCriterion)
                                % visualize the trial
                                figure; hold on;
                                plot(1:2000,data.pupilTracePerTrial(iTrial,:),'-b','LineWidth',1.5);
                                % plot max and min windows
                                plot(maxWindow,data.pupilTracePerTrial(iTrial,maxWindow),'-y','LineWidth',2);
                                plot(minWindow,data.pupilTracePerTrial(iTrial,minWindow),'-y','LineWidth',2);
                                % plot locations in trace that meet velocity criterion
                                exceedsCriterion = data.pupilTracePerTrial(iTrial,:);
                                exceedsCriterion(abs(velocity) <= velocityCriterion) = nan;
                                plot(1:2000,exceedsCriterion,'-r','LineWidth',2);
                                % labels
                                title(['Pupil Size Experiment: ' subjects{iSubj} ' trial#' num2str(iTrial)]);
                                xlabel('Time since stimulus onset (ms)');
                                ylabel('Percent change in pupil Size from baseline');
                                % manual input of decision to include or not
                                disp('Look at the current trial and keep in mind whether the measurement noise will impact the amplitude signal.')
                                decision = input('Should this trial be included in the analysis of pupil amplitude variablity? ');
                                % close figures that are included, otherwise
                                % send to back
                                if decision
                                    % close figures that are included, such that
                                    % only excluded trials remain open in figure
                                    % windows
                                    close gcf
                                else
                                    % send figure to back so that it is not
                                    % confused for a new figure
                                    uistack(gcf,'bottom')
                                end
                                excludeTrial(iSubj,iTrial) = ~logical(decision);
                            end
                        end
                        if ~excludeTrial(iSubj,iTrial)
                            % document amplitude and the time of peak and
                            % trough pupil size
                            amp(iSubj,iTrial) = nanmedian(data.pupilTracePerTrial(iTrial,maxWindow)) - nanmedian(data.pupilTracePerTrial(iTrial,minWindow));
                            timeMax(iSubj,iTrial) = maxIdx;
                            timeMin(iSubj,iTrial) = minIdx;
                        end
                    end
                    
%                     % MPS idea for amplitude (peak median - trough median),
%                     % where peak is at 217 ms for the mean trace and trough
%                     % is at 803 ms (based on overall average trail trace
%                     for iTrial = 1:length(data.pupilTracePerTrial(:,1))
%                         % test: check that this approach is reasonable
%                         figure; hold on
%                         plot(1:2000,data.pupilTracePerTrial(iTrial,:),'-b','LineWidth',1.5);
%                         plot(192:242,data.pupilTracePerTrial(iTrial,192:242),'-y');
%                         plot(778:828,data.pupilTracePerTrial(iTrial,778:828),'-y');
%                         % end test
%                         amp(iSubj,iTrial) = nanmedian(data.pupilTracePerTrial(iTrial,192:242)) - nanmedian(data.pupilTracePerTrial(iTrial,778:828));
%                     end
%                     % put all trials in single matrix in order to get
%                     overall average trial trace
%                     allTrials = [allTrials; data.pupilTracePerTrial];
                    
%                     % (max - min) for each trial
%                     for iTrial = 1:length(data.pupilTracePerTrial(:,1))
%                         amp(iSubj,iTrial) = max(data.pupilTracePerTrial(iTrial,:)) - min(data.pupilTracePerTrial(iTrial,:));
%                     end
                end
            end
        end
    end
    % standard deviation
    ampSD(iSubj) = nanstd(amp(iSubj,:));
    timeMaxSD(iSubj) = nanstd(timeMax(iSubj,:));
    timeMinSD(iSubj) = nanstd(timeMin(iSubj,:));
    % trimmed standard deviation
%     trimmedAmp(iSubj,:) = amp(iSubj,:);
%     trimmedAmp(iSubj,AK_trimIdx(trimmedAmp(iSubj,:),20,'high')) = nan; % trim 20% of the highest amplitudes
%     trimmedAmpSD(iSubj) = nanstd(trimmedAmp(iSubj,:));
end

% % plot amplitudes to check the impact of trimming 20% of highest amplitudes
% figure(1);
% AK_plotSpread_mouseover({ampSD,trimmedAmpSD},{subjects,subjects});

%% saving things

cd('C:\Users\Alex Kale\Documents\MATLAB\CustomMurrayLabScripts\pupilSizeExperiment\amplitudeAnalysisData');

% % figures
% h = get(0,'children');
% for i=1:length(h)
%   saveas(h(i), ['figure10' num2str(i)], 'jpg');
% end

% data
% save('pupilAmpData.mat','subjects','excludeTrial','amp','ampSD','trimmedAmp','trimmedAmpSD');
save('pupilAmpData2.mat','subjects','excludeTrial','amp','ampSD','timeMax','timeMaxSD','timeMin','timeMinSD');

%% plot amplitudes by group to see if there is something there

% index: 1 == ASD; 0 == NT
groupIdx = ~cellfun(@isempty,cellfun(@(x) regexp(x,regexptranslate('wildcard','*G1*')),subjects,'UniformOutput',false));
% plot ampSD
figure;
AK_plotSpread_mouseover({ampSD(groupIdx),ampSD(~groupIdx)},{subjects(groupIdx),subjects(~groupIdx)},'xNames',{'ASD','NT'},'showMM',4);
% plot timeMaxSD
figure;
AK_plotSpread_mouseover({timeMaxSD(groupIdx),timeMaxSD(~groupIdx)},{subjects(groupIdx),subjects(~groupIdx)},'xNames',{'ASD','NT'},'showMM',4);
% plot timeMinSD
figure;
AK_plotSpread_mouseover({timeMinSD(groupIdx),timeMinSD(~groupIdx)},{subjects(groupIdx),subjects(~groupIdx)},'xNames',{'ASD','NT'},'showMM',4);


