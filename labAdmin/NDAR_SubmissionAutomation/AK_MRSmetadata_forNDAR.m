function AK_MRSmetadata_forNDAR( metadata_filename, subjInfo_filename )
%AK_MRSmetadata_forNDAR prepares a .csv file containing appropriately 
% formatted meta-data for all necessary MRS data files and saves this file 
% to the NDAR submission folder with the fullfile name specified in the 
% first function input (metadata_filename). It does this for the subjects 
% whose information is contained in the file designated in the second 
% function input (subjInfo_filename). See AK_prepareNDARsubmission.m 
% documentation for more detailed instructions on the NDAR submission 
% process and file formatting.
%
% HISTORY:
% used to be Spectroscopy_Data_for_NDAR
% Murray Lab 2016
% Created by Alex Kale 6/15/16 

%% establish directory and grab subject info

top_dir = 'L:\MurrayLab\ASD\Data'; % set up base directory
mrs_dirs = {'MRS','MRS2','ftap'}; % MRS folders
mrs_loc_dirs = {'MRS-MT-Left','MRS-MT-Right','MRS-Pariet','MRS-Occip','MRS-Occip2','MRS-Temp'}; % MRS subfolders by location
% for filling in comment field
CommentStr = {'Left MT','Right MT','Parietal','Occipital','Occipital','Temporal'}; % comment string roots

% subjInfo_filename = fullfile(top_dir,'NDARnotes','NDAR_SubjectInfo.xlsx'); % create full filename for NDAR subject info spreadsheet
[~,~,subjInfo] = xlsread(subjInfo_filename); % load subject info
subj_dirs = subjInfo(2:end,1); % designate subject ids

dateFormat = 'mm/dd/yyyy'; % string for date formatting

%% create cell array for data output

NDARdata = cell(6*length(subj_dirs)+1,13); % rows: 10 possible MRS scans per subject (if we repeated the MRS session) * #subjects
NDARdata(1,1:2) = {'cu_spectroscopy','1'};
NDARdata(2,:) = {'subjectkey','src_subject_id','interview_date','interview_age','mrs_raw_data_localizer','mrs_raw_data_ovs','mrs_raw_data_pfile','gender','mrs_raw_data_sdat','mrs_raw_data_spar','mrs_raw_water_sdat','mrs_raw_water_spar','comment'};

%% load data and organize in table, per file, per subject

NDARrow = 3; % designate variable for indexing NDARdata rows
actStr = '*act*';
refStr = '*ref*';

for iS = 1:length(subj_dirs); % cycle through subjects
    subj_dir = fullfile(top_dir,subj_dirs{iS});
    if exist(subj_dir,'dir')==0
        disp(['could not find subject folder for' subj_dirs{iS}]) % message
    elseif exist(subj_dir,'dir')==7
        clear mrsCount mrsIndex
        mrsCount = zeros(length(mrs_loc_dirs),1); % counter for number of each kind of MRS scan per participant
        mrsIndex = cell(length(mrs_loc_dirs),1); % index for NDARdata row for each kind of MRS scan per participant
        for iMRS = 1:length(mrs_dirs); % cycle through potential mrs folders
            for iMRSloc = 1:length(mrs_loc_dirs); % cycle through potential mrs subfolders
                use_dir = fullfile(subj_dir,mrs_dirs{iMRS},mrs_loc_dirs{iMRSloc}); % generate directory in use
                if exist(use_dir,'dir')==0 % mrs directory doesn't exist
                    disp([use_dir ' is not an existing directory']); % message
                else 
                    SDAT_list = dir(fullfile(use_dir,'*.SDAT')); % generate structure containing information about SDAT files
                    SPAR_list = dir(fullfile(use_dir,'*.SPAR')); % generate structure containing information about SPAR files
                    if ~isempty(SDAT_list)
                        NDARdata{NDARrow,1} = subjInfo{iS+1,2}; % store GUID
                        NDARdata{NDARrow,2} = subj_dirs{iS}; % store subject id
                        NDARdata{NDARrow,8} = subjInfo{iS+1,4}; % store gender
                        
                        for iSDAT = 1:length(SDAT_list); % cycle through SDAT files
                           if ~isempty(regexp(SDAT_list(iSDAT).name,regexptranslate('wildcard',actStr),'once'))
                               NDARdata{NDARrow,3} = datestr(SDAT_list(iSDAT).date,dateFormat); % store interview date
                               NDARdata{NDARrow,4} = AK_calculateInterviewAge_forNDAR(subjInfo{iS+1,3}, datestr(SDAT_list(iSDAT).date,dateFormat)); % store age
                               NDARdata{NDARrow,9} = SDAT_list(iSDAT).name; % store filename
                           elseif ~isempty(regexp(SDAT_list(iSDAT).name,regexptranslate('wildcard',refStr),'once'))
                               NDARdata{NDARrow,3} = datestr(SDAT_list(iSDAT).date,dateFormat); % store interview date
                               NDARdata{NDARrow,4} = AK_calculateInterviewAge_forNDAR(subjInfo{iS+1,3}, datestr(SDAT_list(iSDAT).date,dateFormat)); % store age
                               NDARdata{NDARrow,11} = SDAT_list(iSDAT).name; % store filename
                           end
                        end
                    else
                        disp(['SDAT_list is empty for ' use_dir]) % message
                    end
                    if ~isempty(SPAR_list)
                        NDARdata{NDARrow,1} = subjInfo{iS+1,2}; % store GUID
                        NDARdata{NDARrow,2} = subj_dirs{iS}; % store subject id
                        NDARdata{NDARrow,8} = subjInfo{iS+1,4}; % store gender
                        
                        for iSPAR = 1:length(SPAR_list); % cycle through SPAR files
                           if ~isempty(regexp(SPAR_list(iSPAR).name,regexptranslate('wildcard',actStr),'once'))
                               NDARdata{NDARrow,3} = datestr(SPAR_list(iSPAR).date,dateFormat); % store interview date
                               NDARdata{NDARrow,4} = AK_calculateInterviewAge_forNDAR(subjInfo{iS+1,3}, datestr(SPAR_list(iSPAR).date,dateFormat)); % store age
                               NDARdata{NDARrow,10} = SPAR_list(iSPAR).name; % store filename
                           elseif ~isempty(regexp(SPAR_list(iSPAR).name,regexptranslate('wildcard',refStr),'once'))
                               NDARdata{NDARrow,3} = datestr(SPAR_list(iSPAR).date,dateFormat); % store interview date
                               NDARdata{NDARrow,4} = AK_calculateInterviewAge_forNDAR(subjInfo{iS+1,3}, datestr(SPAR_list(iSPAR).date,dateFormat)); % store age
                               NDARdata{NDARrow,12} = SPAR_list(iSPAR).name; % store filename
                           end  
                        end
                    else
                        disp(['SPAR_list is empty for ' use_dir]) % message
                    end
                    mrsCount(iMRSloc) = mrsCount(iMRSloc)+1; % count number of each type of MRS scan
                    mrsIndex{iMRSloc}(mrsCount(iMRSloc)) = NDARrow; % keep track of NDARdata row for each instance of each kind of MRS scan
                    
                    NDARrow = NDARrow+1; % counter
                end
            end
        end
        % fill in comment fields for participant
        for iMRStype = 1:length(mrs_loc_dirs); % cycle through mrs_loc_dirs
            for iComment = 1:mrsCount(iMRStype); % cycle through number of comments needed for each type of MRS scan
                if mrsCount(iMRStype) > 1; % if there is more than one MRS scan of iMRStype
                    NDARdata{mrsIndex{iMRStype}(iComment),13} = [CommentStr{iMRStype} ' ' num2str(iComment)]; % store comment string with iComment appended
                else
                    NDARdata(mrsIndex{iMRStype}(iComment),13) = CommentStr(iMRStype); % store comment string
                end
            end
        end
    end
 end

%% find and delete blank rows

isEmptyIndex = zeros(length(NDARdata(:,1)),1);

for iNDAR = 1:length(NDARdata(:,1)); % cycle through rows of NDARdata
    if mean(cellfun(@isempty,NDARdata(iNDAR,:))) == 1; % is the entire row empty?
        isEmptyIndex(iNDAR) = 1;
    end
end

isEmptyIndex = logical(isEmptyIndex); % make logical

NDARdata(isEmptyIndex,:) = []; % clear empty rows

%% write to csv file

cell2csv(metadata_filename,NDARdata);

end
