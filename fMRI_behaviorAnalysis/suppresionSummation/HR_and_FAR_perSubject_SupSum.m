% HR_and_FAR_perSubject_SupSum
% MurrayLab_2016
% Created by AMK 11/2/16

%% directories

top_dir = 'L:\MurrayLab\ASD\SuppressionSummation';
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab\SuppressionSummation';
subjects = {'S_AMK_20160401','S_AVF_20160427','S_AVF_20160510','S_AW_20160420','S_BK_20160517','S_DP_20160427','S_JCM_20160426','S_MN_20160523','S_MPS_20160325','S_RM_20160414','S_SOM_20160414'};

%% create matrix containing data

avgHR_FAR = nan(length(subjects),2);

for iS = 1:length(subjects)
    % create directory and look for it
    use_dir = fullfile(top_dir,subjects{iS},'log_files');
    if exist(use_dir,'dir')
        % find logfiles
        log_list = dir(fullfile(use_dir,'*.log'));
        % sort by number following 'initials_'
        underscoreIdx = strfind(subjects{iS},'_');
        log_list = AK_sortStruct(log_list,1,subjects{iS}(3:underscoreIdx(2))); % initials always start at 3rd position and end at the second underscore
        % setup
        clear HR FAR
        HR = nan(size(log_list));
        FAR = nan(size(log_list));
        for iL = 1:length(log_list)
            % parse logfile
            clear a
            [a,~,~,~] = AK_SupSum_logfileParse(fullfile(use_dir,log_list(iL).name));
            % record hit rate and false alarm rate per logfile
            HR(iL) = a.hitRate;
            FAR(iL) = a.falseAlarmRate;
        end
        % calculate average hit rate and false alarm rate
        avgHR_FAR(iS,1) = nanmean(HR);
        avgHR_FAR(iS,2) = nanmean(FAR);
    else
        avgHR_FAR(iS,1) = nan;
        avgHR_FAR(iS,2) = nan;
    end
end

%% save

save(fullfile(home_dir,'avgHRandFAR.mat'),'avgHR_FAR','subjects');
