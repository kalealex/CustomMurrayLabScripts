% ftap Log File Processing Pipeline
% Murray Lab_2015
% Created by Alex Kale 11/30/15
% Some code from Michael-Paul

%% establish directory by subjects and fMRI sessions

top_dir = 'L:\MurrayLab\ASD\Data'; % set up base directory
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab'; % set directory for saving data and figures
subj_dir = {'G313'}; % designate folders for invidual participants 
fmri_dir = {'ftap'}; % designate subfolders
conds = {'standard', 'variable'}; % create list of conditions for figure creation

%% parse log files and store in data structure

for iS = 1:length(subj_dir); % cycle through subjects
    for iF = 1:length(fmri_dir); % cycle through fmri sessions
        use_dir = fullfile(top_dir,subj_dir{iS},fmri_dir{iF},'log_files'); % generate directory in use
        log_list = dir(fullfile(use_dir,'*.log')); % generate structure containing information about log files
        if ~isempty(log_list)
            for iL = 1:length(log_list) % cycle through log files
                clear accuracy resptimes timepts
                mat_name = [fullfile(use_dir,log_list(iL).name(1:end-3)) 'mat']; % create string to name .mat file which saves timepts array
                    if ~exist(mat_name,'file') % check whether or not function output has already been saved
                        disp(['parsing ' fullfile(use_dir,log_list(iL).name)]) % message
                        [accuracy,resptimes,timepts] = AK_GABAinASD_ftap_logfileParse(fullfile(use_dir,log_list(iL).name)); % run log file through parsing function
                        save(mat_name,'accuracy','resptimes','timepts'); % save function output
                    else
                        disp(['loading ' mat_name]) % message
                        load(mat_name,'accuracy','resptimes','timepts') % load file containing timepts array
                    end    
%                 data.logFile{iS,iF,iL} = log_list(iL).name; % store name of logfile
%                 data.hitRateStd(iS,iF,iL) = accuracy(1); % store hit rate information for standard stim in data structure
%                 data.hitRateVar(iS,iF,iL) = accuracy(2); % store hit rate information for variable stim in data structure
%                 data.missRateStd(iS,iF,iL) = accuracy(3); % store miss rate information for standard stim in data structure
%                 data.missRateVar(iS,iF,iL) = accuracy(4); % store miss rate information for variable stim in data structure
%                 data.errorRate(iS,iF,iL) = accuracy(5); % store error rate information in data structure
%                 data.respCountAll(iS,iF,iL) = accuracy(6);
%                 data.resptimesAllM(iS,iF,iL) = nanmean(resptimes.all);
%                 data.resptimesAllSD(iS,iF,iL) = nanstd(resptimes.all);
%                 data.respCountAllStd(iS,iF,iL) = length(resptimes.allStd);
%                 data.resptimesAllStdM(iS,iF,iL) = nanmean(resptimes.allStd);
%                 data.resptimesAllStdSD(iS,iF,iL) = nanstd(resptimes.allStd);
%                 data.respCountAllVar(iS,iF,iL) = length(resptimes.allVar);
%                 data.resptimesAllVarM(iS,iF,iL) = nanmean(resptimes.allVar);
%                 data.resptimesAllVarSD(iS,iF,iL) = nanstd(resptimes.allVar);
%                 data.resptimesAll(iS,iF,iL,1:length(resptimes.all)) = resptimes.all; % store response times in fourth dimension of matrix
%                 data.resptimesAllStd(iS,iF,iL,1:length(resptimes.allStd)) = resptimes.allStd;
%                 data.resptimesAllVar(iS,iF,iL,1:length(resptimes.allVar)) = resptimes.allVar;
%                 data.resptimesHitsStd(iS,iF,iL,1:length(resptimes.hitsStd)) = resptimes.hitsStd;
%                 data.resptimesHitsVar(iS,iF,iL,1:length(resptimes.hitsVar)) = resptimes.hitsVar;
%                 data.timepts{iS,iF,iL} = timepts; % store timepts array
            end
        end
    end
end
% save data structure
% data_dir = fullfile(home_dir,'Behavioral_Data',subj_dir{iS},fmri_dir{iF}); % set directory for saving data
% data_name = fullfile(data_dir,'ftap_Behavioral_Data.mat'); % create string to name .mat file which saves data
% if ~exist(data_name,'file') % check whether or data have already been saved
%     disp(['saving ' data_name]) % message
%     save(data_name,'data'); % save data
% end

%% figures 

% for s = 1:length(subj_dir); % cycle through subjects
%     for f = 1:length(fmri_dir); % cycle through fmri sessions
%         use_dir = fullfile(top_dir,subj_dir{iS},fmri_dir{f},'log_files'); % generate directory in use
%         log_list = dir(fullfile(use_dir,'*.log')); % generate structure containing information about log files
%         for c = 1:length(conds);
%             fig_name = ['Response Times per Scan for ' fmri_dir{f} ', ' subj_dir{s} ': ' conds{c}] ;
%             figure('Name',fig_name,'NumberTitle','off') % create resptimes histograms for standard stim, for each subject
%             if ~isempty(log_list)
%                 for l = 1:length(log_list) % cycle through log files
%                     % set subplot dimensions
%                     if length(log_list) > 7
%                         subplot(2,4+mod(length(log_list),4),l);
%                     elseif length(log_list) > 5
%                         subplot(2,3+mod(length(log_list),3),l);
%                     else
%                         subplot(2,2+mod(length(log_list),2),l);
%                     end
%                     % make histogram for either standard or variable stim
%                     if c ==1;
%                         hist(squeeze(data.resptimesAllStd(s,f,l,:)));
%                     elseif c==2;
%                         hist(squeeze(data.resptimesAllVar(s,f,l,:)));
%                     end
%                     title(log_list(l).name(1:end-4));
%                     xlabel('Response Times (1/10000 sec)');
%                     ylabel('Frequency'); 
%     %                 set(gca,'XTick',0:5000:max(squeeze(data.resptimesAll(s,f,l,:))));
% 
%     %                 descriptives.hitRatemean(s) = mean(reshape(data.hitRate(s,:,:),[length(fmri_dirs)*length(data.hitRate(s,:,:)),1])); % create structure of descriptive statistics per subject
%     %                 descriptives.hitRateSD(s) = std(reshape(data.hitRate(s,:,:),[length(fmri_dirs)*length(data.hitRate(s,:,:)),1]));
%     %                 descriptives.resptimesAllmean(s) = nanmean(reshape(data.resptimesAll(s,:,:,:),[length(fmri_dirs)*50*100,1]));
%     %                 descriptives.resptimesAllSD(s) = nanstd(reshape(data.resptimesAll(s,:,:,:),[length(fmri_dirs)*50*100,1]));
%                 end
%             end
%             % save figure
%             fig_dir = fullfile(home_dir,'Behavioral_Data',subj_dir{s},fmri_dir{f}); % set directory for saving current figure
%             fig_name = fullfile(fig_dir,[conds{c} '_ftapResponseTimes.fig']); % create string to name .fig file which saves current figure
%             if ~exist(fig_name,'file') % check whether or not function output has already been saved
%                 disp(['saving ' fig_name]) % message
%                 saveas(gcf,fig_name,'fig'); % save figure
%             end
%         end
%     end
% end

% descriptives.hitRateGrandmean = mean(descriptives.hitRatemean(1:length(subj_dirs))); % Calculate Grand descriptives (across subjects
% descriptives.hitRateGrandSD = std(reshape(data.hitRate(:,:,:),[length(subj_dirs)*length(fmri_dirs)*length(data.hitRate(:,:,:)),1]));
% descriptives.resptimesAllGrandmean = mean(descriptives.resptimesAllmean(1:length(subj_dirs)));
% descriptives.resptimesAllGrandSD = nanstd(reshape(data.resptimesAll(:,:,:,:),[length(subj_dirs)*length(fmri_dirs)*50*100,1]));

