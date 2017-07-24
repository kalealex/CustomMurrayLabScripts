% Fill_GABAinASD_SRS_RMET_SENSORY_ADOS_DataStructure_fromXLSfile
% MurrayLab 2017
% created by AMK 170628

% load data from speadsheet:
% this spreadsheet might need to be updated by copying in recent entries
% from REDCap Report: Michael-Paul, ADOS, SRS, RMET, & Sensory Profile, 20170113
[~,~,data] = xlsread('L:\MurrayLab\ASD\Data\ADOS, SRS, RMET, & Sensory Profile.xlsx');

% make each column a vector field of the structure symtom_scores
symptom_scores.subj_num = cellfun(@str2double,cellfun(@(x) AK_eraseSubstring(x,'G'),data(4:end,1),'UniformOutput',false));
symptom_scores.SRS_total = cell2mat(data(4:end,2));
symptom_scores.RMET_total = cell2mat(data(4:end,3));
symptom_scores.sensory_low_registration = cell2mat(data(4:end,4));
symptom_scores.sensory_sensation_seeking = cell2mat(data(4:end,5));
symptom_scores.sensory_sensitivity = cell2mat(data(4:end,6));
symptom_scores.sensory_avoiding = cell2mat(data(4:end,7));
symptom_scores.ADOS_social_affect = cell2mat(data(4:end,8));
symptom_scores.ADOS_restricted_repetative_behavior = cell2mat(data(4:end,9));
symptom_scores.ADOS_total = cell2mat(data(4:end,10));

% save
% save(['L:\MurrayLab\ASD\Data\Analysis_scripts\GABAinASD_SymptomScores' datestr(now,'yyyymmdd') '.mat'],'symptom_scores');