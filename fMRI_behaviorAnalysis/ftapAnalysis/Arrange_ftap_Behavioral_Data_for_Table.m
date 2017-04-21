% Arrange_Behavioral_Data_for_Table
% Murray Lab_2015
% Created by Alex Kale 11/30/15

%% establish directory to retrieve data struct by subjects and fMRI sessions

home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab\Behavioral_Data'; % set directory for saving figures
subj_dir = 'G313'; % designate folder for invidual participants 
fmri_dir = 'ftap'; % designate ftap folder

cd(home_dir); % set directory

%% load data structure

filename = fullfile(home_dir,subj_dir,fmri_dir,'ftap_Behavioral_Data.mat');
if exist(filename,'file') % attempt to load data struct by both possible names
    load(filename,'data');
else
    disp('file does not exist');
end

%% establish ordered list of fields and store info from data struct in table

fields = fieldnames(data); % create cell array of field names in order

table = cell(length(data.logFile(1,1,:))+1,length(fields)-6); % set up a cell array the proper size to store all the data
table(1,1:end) = fields(1:end-6); % set table column labels

for f = 1:length(fields)-6; % cycle through fields 'fMRIsession' - 'resptimesAllSD'
    row = 2; % create/set counter for table row indexing
    if ~isempty(data.logFile(1,1,:));
        for l = 1:length(data.logFile(1,1,:)); % cycle through list of .log filenames
            if f == 1
                table(row,f) = data.logFile(1,1,l); % store log file name
                row = row+1; % add one to table row counter
            else
                table{row,f} = getfield(data,fields{f},{1,1,l}); % store data for  current .log file, for current fmri session, for current field
                row = row+1; % add one to table row counter
            end
        end
    end
end

posNoFile = cellfun(@isempty,table(:,1)); % find rows in which there is no log file
table(posNoFile,:) = []; % delete rows where there is no log file

%% save table of data

table_dir = fullfile(home_dir,subj_dir,fmri_dir); % set directory for saving table
table_name = [fullfile(table_dir,subj_dir) '_ftap_Behavioral_Data_table.mat']; % create string to name .mat file which saves table
if ~exist(table_name,'file') % check whether or table has already been saved
    disp(['saving ' table_name]) % message
    save(table_name,'table'); % save table
end



