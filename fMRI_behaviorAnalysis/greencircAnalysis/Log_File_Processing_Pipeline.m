% Log File Processing Pipeline
% Murray Lab_2015
% Created by Alex Kale 9/23/15
% Some code from Michael-Paul

%% establish directory by subjects and fMRI sessions

top_dir = 'L:\MurrayLab\ASD\Data'; % set up base directory
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab'; % set directory for saving data and figures
subj_dirs = {'G316'}; % designate folders for invidual participants 
fmri_dirs = {'fMRI1','fMRI2'}; % designate subfolders

%% parse log files and store in data structure

for iS = 1:length(subj_dirs); % cycle through subjects
    for iF = 1:length(fmri_dirs); % cycle through fmri sessions
        use_dir = fullfile(top_dir,subj_dirs{iS},fmri_dirs{iF},'log_files'); % generate directory in use
        log_list = dir(fullfile(use_dir,'*.log')); % generate structure containing information about log files
        if ~isempty(log_list)
            for iL = 1:length(log_list) % cycle through log files
                clear accuracy resptimes timepts
                mat_name = [fullfile(use_dir,log_list(iL).name(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
%                     if ~exist(mat_name,'file') % check whether or not function output has already been saved
                        disp(['parsing ' fullfile(use_dir,log_list(iL).name)]) % message
                        [accuracy,resptimes,timepts] = AK_GABAinASD_logfileParse(fullfile(use_dir,log_list(iL).name)); % run log file through parsing function
                        save(mat_name,'accuracy','resptimes','timepts'); % save function output
%                     else
%                         disp(['loading ' mat_name]) % message
%                         load(mat_name,'accuracy','resptimes','timepts') % load file containing timepts array
%                     end
                data.fMRIsession{iS,iF} = fmri_dirs(iF);    
                data.logFile{iS,iF,iL} = log_list(iL).name; % store name of logfile
                data.hitRate(iS,iF,iL) = accuracy(1); % store hit rate information in data structure
                data.missRate(iS,iF,iL) = accuracy(2); % store miss rate information in data structure
                data.errorRate(iS,iF,iL) = accuracy(3); % store error rate information in data structure
                data.respCountAll(iS,iF,iL) = accuracy(4);
                data.presentationCountAll(iS,iF,iL) = accuracy(5);
                data.resptimesAllM(iS,iF,iL) = nanmean(resptimes.all);
                data.resptimesAllSD(iS,iF,iL) = nanstd(resptimes.all);
                data.resptimesAll(iS,iF,iL,1:length(resptimes.all)) = resptimes.all; % store response times in fourth dimension of matrix
                data.resptimesHits(iS,iF,iL,1:length(resptimes.hits)) = resptimes.hits;
                data.timepts{iS,iF,iL} = timepts; % store timepts array
            end
        end
    end
end
% save data structure
data_dir = fullfile(home_dir,'Behavioral_Data',subj_dirs{iS}); % set directory for saving data
data_name = fullfile(data_dir,'fMRI_Behavioral_Data.mat'); % create string to name .mat file which saves data
% if ~exist(data_name,'file') % check whether or data have already been saved
    disp(['saving ' data_name]) % message
    save(data_name,'data'); % save data
% end

%% figures 

for s = 1:length(subj_dirs); % cycle through subjects
    for f = 1:length(fmri_dirs); % cycle through fmri sessions
        use_dir = fullfile(top_dir,subj_dirs{iS},fmri_dirs{f},'log_files'); % generate directory in use
        log_list = dir(fullfile(use_dir,'*.log')); % generate structure containing information about log files
        fig_name = ['Response Times per Scan for ' fmri_dirs{f} ', ' subj_dirs{iS}] ;
        figure('Name',fig_name,'NumberTitle','off') % create resptimes histograms for each subject
        if ~isempty(log_list)
            for l = 1:length(log_list) % cycle through log files
                if length(log_list) > 7
                    subplot(2,4+mod(length(log_list),4),l);
                elseif length(log_list) > 5
                    subplot(2,3+mod(length(log_list),3),l);
                else
                    subplot(2,2+mod(length(log_list),2),l);
                end
                hist(squeeze(data.resptimesAll(s,f,l,:)));
                title(log_list(l).name(1:end-4));
                xlabel('Response Times (1/10000 sec)');
                ylabel('Frequency'); 
%                 set(gca,'XTick',0:5000:max(squeeze(data.resptimesAll(s,f,l,:))));
                
%                 descriptives.hitRatemean(s) = mean(reshape(data.hitRate(s,:,:),[length(fmri_dirs)*length(data.hitRate(s,:,:)),1])); % create structure of descriptive statistics per subject
%                 descriptives.hitRateSD(s) = std(reshape(data.hitRate(s,:,:),[length(fmri_dirs)*length(data.hitRate(s,:,:)),1]));
%                 descriptives.resptimesAllmean(s) = nanmean(reshape(data.resptimesAll(s,:,:,:),[length(fmri_dirs)*50*100,1]));
%                 descriptives.resptimesAllSD(s) = nanstd(reshape(data.resptimesAll(s,:,:,:),[length(fmri_dirs)*50*100,1]));
            end
        end
%         % save figure
%         fig_dir = fullfile(home_dir,'Behavioral_Data',subj_dirs{s},fmri_dirs{f}); % set directory for saving current figure
%         fig_name = fullfile(fig_dir,'ResponseTimes.fig'); % create string to name .fig file which saves current figure
%         if ~exist(fig_name,'file') % check whether or not function output has already been saved
%             disp(['saving ' fig_name]) % message
%             saveas(gcf,fig_name,'fig'); % save figure
%         end
    end
end

% descriptives.hitRateGrandmean = mean(descriptives.hitRatemean(1:length(subj_dirs))); % Calculate Grand descriptives (across subjects
% descriptives.hitRateGrandSD = std(reshape(data.hitRate(:,:,:),[length(subj_dirs)*length(fmri_dirs)*length(data.hitRate(:,:,:)),1]));
% descriptives.resptimesAllGrandmean = mean(descriptives.resptimesAllmean(1:length(subj_dirs)));
% descriptives.resptimesAllGrandSD = nanstd(reshape(data.resptimesAll(:,:,:,:),[length(subj_dirs)*length(fmri_dirs)*50*100,1]));

