% LZ_FixationDifferences_DrugvsPlacebo
% MurrayLab 2017
% created by AMK 6/2/17

%% load data

cd('L:\MurrayLab\Lorazepam\analysis_code');
load('LorazepamMotion_EyeTrackingData_driftcorrDistM_20170602.mat'); % output from the GUI :D

conditions = key(1).experiment.conditions;

%% extract statistics

% preallocate
drugFixation = cell(0);
placeboFixation = cell(0);
pRankSum = nan(size(conditions));
% extract a set of average drift corrected distances from fixation for each
% condition (subjects * runs)
for iCond = 1:length(conditions)
    drugFixation{iCond} = squeeze(exportData(1,iCond,1,:,1,:));
    placeboFixation{iCond} = squeeze(exportData(1,iCond,1,:,2,:));
    pRankSum(iCond) = signrank(nanmean(drugFixation{iCond},2) - nanmean(placeboFixation{iCond},2));
end

%% big picture stats

overallPvalue = pRankSum(8); % across all conditions

avgOverallDifferencePerSubj = nanmean(drugFixation{8},2) - nanmean(placeboFixation{8},2); % degrees visual angle; averaged across runs within each subject prior to taking difference

avgDifferenceAcrossSubjects = mean(avgOverallDifferencePerSubj); % degrees visual angle

effectSize = avgDifferenceAcrossSubjects/std(avgOverallDifferencePerSubj);

%% impact on thresholds

% get threshold data
load('LZdata_20170113.mat');

% preallocate
drugThresholds = cell(0);
placeboThresholds = cell(0);
% extract thresholds into similarly formatted cell array
for iCond = 1:length(conditions)-1 % thresholds are not aggregated across trials the way that eyetracking statistics are, so travers one fewer conditions
    drugThresholds{iCond} = squeeze(LZssMotion.thresholds(:,1,iCond,:));
    placeboThresholds{iCond} = squeeze(LZssMotion.thresholds(:,2,iCond,:));
end

% correlation analysis:

% preallocate
fixDiff = cell(0);
thresDiff = cell(0);
rCorr = cell(0);
pCorr = cell(0);
% process conditions separately
for iCond = 1:length(conditions)-2 % exclude catch trials
    % find differences between sessions, averaging across runs
    fixDiff{iCond} = nanmean(drugFixation{iCond},2) - nanmean(placeboFixation{iCond},2);
    thresDiff{iCond} = nanmean(drugThresholds{iCond},2) - nanmean(placeboThresholds{iCond},2);
    % correlation
    [rCorr{iCond},pCorr{iCond}] = corrcoef(fixDiff{iCond}, thresDiff{iCond});
end
