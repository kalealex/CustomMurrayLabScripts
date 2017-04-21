% GABAinASD_fMRI_Noise_CRF_MT
% MurrayLab_2016
% created by AMK 10/14/16

% switch:
AV = 0; % averaged across trials or not (boolean)

%% directory components

% root dirs
top_dir = 'L:\MurrayLab\ASD\Data';
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab\GABAinASD_fMRI_Noise';
% all subjects; won't work for subjects with unanalyzed fMRI data
subj_dirs = {'G101','G102','G103','G104','G105','G106','G107','G109','G110','G111','G112','G114','G117','G119','G307','G310','G311','G312','G313','G314','G315','G316','G317','G318','G319','G320','G322','G324','G326','G327','G328','G329','G330','G332','G336','G338','G345'};
% session subfolders
fmri_dirs = {'fMRI1','fMRI2'};
% mat files to load
mat_file_MTleft = 'CRF_MT_MT*left.mat';
mat_file_MTright= 'CRF_MT_MT*right.mat';

% key for columns of noiseData matrix
noiseKey = {'subj#','LC_LH_M','LC_LH_ST','LC_RH_M','LC_RH_ST','HC_LH_M','HC_LH_ST','HC_RH_M','HC_RH_ST'};

% disqualified subjects
dqSubjects = {'G103','G117','G310','G326'};

%% query directories to fill in noiseData matrix

% preallocate
noiseData = nan(length(subj_dirs)*length(fmri_dirs),length(noiseKey));
% set counter
ndRow = 1;

for iS = 1:length(subj_dirs)
    for iF = 1:length(fmri_dirs)
        % create directory in use and check that it exists
        clear use_dir
        use_dir = fullfile(top_dir,subj_dirs{iS},fmri_dirs{iF},'RESULTS','data');
        if exist(use_dir,'dir')
            % check that mat files to load exist:
            % for left MT
            clear mat_file AVLowContrastData LowContrastData AVHighContrastData HighContrastData
            mat_file = dir(fullfile(use_dir,mat_file_MTleft));
            % check for multiple matching filenames
            if length(mat_file)>1 % multiple files
                disp(['*** found multiple mat files for ' subj_dirs{iS} ': ' fmri_dirs{iF}]); % message
                for iM = 1:length(mat_file)
                    disp([num2str(iM) ': ' mat_file(iM).name])
                end
                matIdx = input('Pick a file index number:');
                mat_file = mat_file(matIdx).name;
                mat_file = fullfile(use_dir,mat_file);
            elseif length(mat_file)==1 % one file
                mat_file = mat_file(1).name;
                mat_file = fullfile(use_dir,mat_file);
            else % no file
                mat_file = [];
            end 
            if exist(mat_file,'file')
                % load data
                load(mat_file);
                % calculate noise
                switch AV
                    case 1
                        LC_LH_M = nanmean(AVLowContrastData); % low contrast, left hemi
                        LC_LH_ST = nanstd(AVLowContrastData);
                        HC_LH_M = nanmean(AVHighContrastData); % high contrast, left hemi
                        HC_LH_ST = nanstd(AVHighContrastData);
                    case 0
                        LC_LH_M = mean(nanmean(LowContrastData(:,5:9))); % low contrast, left hemi
                        LC_LH_ST = mean(nanstd(LowContrastData(:,5:9)));
                        HC_LH_M = mean(nanmean(HighContrastData(:,5:9))); % high contrast, left hemi
                        HC_LH_ST = mean(nanstd(HighContrastData(:,5:9)));
                end
%                 % test
%                 figure
%                 plot(mean(LowContrastData),'b-')
%                 hold on;
%                 plot(mean(HighContrastData),'r-')
%                 title([subj_dirs{iS} ': ' fmri_dirs{iF} ' left MT'])
            else
                % missing values = nan
                LC_LH_M = nan; 
                LC_LH_ST = nan;
                HC_LH_M = nan; 
                HC_LH_ST = nan;
            end
            % for right MT
            clear mat_file AVLowContrastData LowContrastData AVHighContrastData HighContrastData
            mat_file = dir(fullfile(use_dir,mat_file_MTright));
            % check for multiple matching filenames
            if length(mat_file)>1 % multiple files
                disp(['*** found multiple mat files for ' subj_dirs{iS} ': ' fmri_dirs{iF}]); % message
                for iM = 1:length(mat_file)
                    disp([num2str(iM) ': ' mat_file(iM).name])
                end
                matIdx = input('Pick a file index number:');
                mat_file = mat_file(matIdx).name;
                mat_file = fullfile(use_dir,mat_file);
            elseif length(mat_file)==1 % one file
                mat_file = mat_file(1).name;
                mat_file = fullfile(use_dir,mat_file);
            else % no file
                mat_file = [];
            end 
            if exist(mat_file,'file')
                % load data
                load(mat_file);
                % calculate noise
                switch AV
                    case 1
                        LC_RH_M = nanmean(AVLowContrastData); % low contrast, right hemi
                        LC_RH_ST = nanstd(AVLowContrastData);
                        HC_RH_M = nanmean(AVHighContrastData); % high contrast, right hemi
                        HC_RH_ST = nanstd(AVHighContrastData);
                    case 0
                        LC_RH_M = mean(nanmean(LowContrastData(:,5:9))); % low contrast, right hemi
                        LC_RH_ST = mean(nanstd(LowContrastData(:,5:9)));
                        HC_RH_M = mean(nanmean(HighContrastData(:,5:9))); % high contrast, right hemi
                        HC_RH_ST = mean(nanstd(HighContrastData(:,5:9)));
                end
%                 % test
%                 figure
%                 plot(mean(LowContrastData),'b-')
%                 hold on;
%                 plot(mean(HighContrastData),'r-')
%                 title([subj_dirs{iS} ': ' fmri_dirs{iF} ' right MT'])
            else
                % missing values = nan
                LC_RH_M = nan; 
                LC_RH_ST = nan;
                HC_RH_M = nan; 
                HC_RH_ST = nan;
            end
        else
           % missing values = nan
            LC_LH_M = nan; 
            LC_LH_ST = nan;
            HC_LH_M = nan; 
            HC_LH_ST = nan; 
            LC_RH_M = nan; 
            LC_RH_ST = nan;
            HC_RH_M = nan; 
            HC_RH_ST = nan;
        end
        % store noise data in matrix
        noiseData(ndRow,:) = [str2double(subj_dirs{iS}(2:end)) LC_LH_M LC_LH_ST LC_RH_M LC_RH_ST HC_LH_M HC_LH_ST HC_RH_M HC_RH_ST];
        clear LC_LH_M LC_LH_ST LC_RH_M LC_RH_ST HC_LH_M HC_LH_ST HC_RH_M HC_RH_ST
        % advance counter
        ndRow = ndRow+1; 
    end
end

%% save 

switch AV
    case 1
        save_file = fullfile(home_dir,'fMRI_Noise_CRF_MT.mat');
    case 0
        save_file = fullfile(home_dir,'fMRI_Noise_CRF_MT_nonAV.mat');
end
save(save_file,'noiseData','noiseKey','dqSubjects');
