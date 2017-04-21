% Analyze_All_fMRI_EyeTracking_Data_GABAinASD
% Murray_Lab_2016
% Created by AMK on 7/12/16, modefies to add block data 7/28/16

%% designate directories 

top_dir = 'L:\MurrayLab\ASD\Data';
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab';

%% load eyetracking quality table to find subjects to include in analysis

load(fullfile(top_dir,'EyetrackingQuality.mat'));
subjects = qualityTable(2:end,1); % full list
subjects(cellfun(@isempty,subjects)) = []; % remove empty cells
subjects = unique(subjects); % remove repetitions from list

% to avoid repeat analysis
finishedSubj = {'G101','G102','G103','G104','G105','G106','G107','G109','G110','G111','G112','G114','G119','G121',...
    'G124','G306','G307','G310','G311','G312','G313','G314','G315','G316','G317','G318','G319','G320','G322','G324',...
    'G326','G327','G328','G329','G330','G331','G332','G336','G338','G340','G342','G345'};
subjects = setdiff(subjects,finishedSubj);

%% run function and save data

for iS = 1:length(subjects) % cycle through subjects 
    % parse eye tracking data and save data structure for each subject
    % individually
    clear eyetracking_data_filename eyetracking_data_dir data
    eyetracking_data_filename = [subjects{iS} '_fMRI_EyeTracking_Data.mat'];
    eyetracking_data_dir = fullfile(top_dir,subjects{iS},eyetracking_data_filename);
%     if ~exist(eyetracking_data_dir,'file') % look for saved data
        data = AK_GABAinASD_fMRI_EvaluateFixation(subjects(iS));
        save(eyetracking_data_dir,'data');
%     else
%         disp(['loading ' eyetracking_data_dir]) % message
%         load(eyetracking_data_dir,'data');
%     end
end
