% Get Demographics and Symptom Scores for GABAinASD
% MurrayLab 2017
% created by AMK 170628, updated to include demographics 20170814

% History: used to be
% 'Fill_GABAinASD_SRS_RMET_SENSORY_ADOS_DataStructure_fromXLSfile.m'

% ATTENTION: For updated data, copy and paste data from REDCap Report "Alex
% Kale, Demographics and Symptom Scores, 20170814" into the file
% 'L:\MurrayLab\ASD\Data\Demographics and Symptom Scores.xlsx' prior to
% running this script. Also, note that the 'demographics' and
% 'symptom_scores' structures must be saved manually, in accord with data
% management conventions for this level of analysis.

% load data from speadsheet:
% this spreadsheet might need to be updated by copying in recent entries
% from REDCap Report: Alex Kale, Demographics and Symptom Scores, 20170814
% used to grab symptom scores from Michael-Paul, ADOS, SRS, RMET, & Sensory Profile, 20170113
[~,~,data] = xlsread('L:\MurrayLab\ASD\Data\Demographics and Symptom Scores.xlsx');

%% demographics

% make the vectors in columns 1:4 into field of the structure demographics
demographics.subj_number = cellfun(@str2double,cellfun(@(x) AK_eraseSubstring(x,'G'),data(4:end,1),'UniformOutput',false));
demographics.interview_age_in_months = cell2mat(data(4:end,2));
demographics.sex = cell2mat(cellfun(@(x) AK_getNumericChars(x),data(4:end,3),'UniformOutput',false));
demographics.non_verbal_IQ = cell2mat(data(4:end,4));
demographics.sex_key = '0 = male, 1 = female';

%% symptom scores

% make the vectors in columns 1 and 5:13 into field of the structure
% symptom_scores
symptom_scores.subj_num = cellfun(@str2double,cellfun(@(x) AK_eraseSubstring(x,'G'),data(4:end,1),'UniformOutput',false));
symptom_scores.SRS_total = cell2mat(data(4:end,5));
symptom_scores.RMET_total = cell2mat(data(4:end,6));
symptom_scores.sensory_low_registration = cell2mat(data(4:end,7));
symptom_scores.sensory_sensation_seeking = cell2mat(data(4:end,8));
symptom_scores.sensory_sensitivity = cell2mat(data(4:end,9));
symptom_scores.sensory_avoiding = cell2mat(data(4:end,10));
symptom_scores.ADOS_social_affect = cell2mat(data(4:end,11));
symptom_scores.ADOS_restricted_repetative_behavior = cell2mat(data(4:end,12));
symptom_scores.ADOS_total = cell2mat(data(4:end,13));

