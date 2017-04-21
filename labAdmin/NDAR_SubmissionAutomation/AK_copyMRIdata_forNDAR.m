function AK_copyMRIdata_forNDAR( destination_dir, subjInfo_filename )
%AK_copyMRIdata_forNDAR copies all necessary fMRI and MRS data files from 
% the L drive to the NDAR submission folder specified in the first function
% input (destination_dir). It does this for the subjects whose information
% is contained in the file designated in the second function input
% (subjInfo_filename). See AK_prepareNDARsubmission.m documentation for
% more detailed instructions on the NDAR submission process and file
% formatting.
% 
% HISTORY:
% used to be fMRI_and_MRS_Copyfile_for_NDAR.m
% Murray Lab 2015
% Created by Alex Kale 1/11/16
% Edited 6/16/16
% Turned into function 12/14/16

%% establish directory by subjects

top_dir = 'L:\MurrayLab\ASD\Data'; % set up base directory
fMRI_dirs = {'fMRI1','fMRI2','ftap'}; % designate folders for fMRI sessions
scan_dirs = {'3danat','set1','set2','set3','set4','set5','set6','set7','set8','set9'}; % designate folders for scans
MRS_dirs = {'MRS','MRS1','MRS2','ftap'}; % designate folders for MRS sessions
subMRS_dirs = {'MRS-MT-Left','MRS-MT-Right','MRS-Pariet','MRS-Occip','MRS-Occip2','MRS-Temp'}; % designate subfolders for MRS
% destination_dir = 'L:\MurrayLab\ASD\Data\NDARnotes\NDARsubmission_June2016'; % designate destination directory
dynStr = {'MPRAGE','125Dyn','80Dyn','170Dyn'};

% subjInfo_filename = fullfile(top_dir,'NDARnotes','NDAR_SubjectInfo.xlsx'); % create full filename for NDAR subject info spreadsheet
[~,~,subjInfo] = xlsread(subjInfo_filename); % load subject info
subj_dirs = subjInfo(2:end,1); % designate subject ids

%% create cell array for data output

NDARdata = cell(10*length(subj_dirs)+1,10);
NDARdata(1,:) = {'subject','fMRI_MRS','scan','date','PARs','RECs','actSPAR','actSDAT','refSPAR','refSDAT'};

%% load data and organize in table, per file, per subject

NDARrow = 2; % designate variable for indexing NDARdata rows
nScans = zeros(length(subj_dirs),1); % create vector for storing number of scans per subject
nMRSs = zeros(length(subj_dirs),1); % create vector for storing number of MRSs per subject

for iS = 1:length(subj_dirs); % cycle through subjects
    subj_dir = fullfile(top_dir,subj_dirs{iS});
    if exist(subj_dir,'dir')==0
        disp(['could not find subject folder for' subj_dirs{iS}]) % message
    elseif exist(subj_dir,'dir')==7
        for iF = 1:length(fMRI_dirs); % cycle through fMRI folders
            clear nScansExpected iScan % set the number of expected scans
            if iF == 3 % for ftap 
                nScansExpected = 4;
            else % for fMRI 1 & 2
                nScansExpected = length(scan_dirs);
            end
            iScan = 1; % reset iScan
            while iScan <= nScansExpected; % cycle through scans
                use_dir = fullfile(subj_dir,fMRI_dirs{iF},scan_dirs{iScan}); % generate directory in use
                if exist(use_dir,'dir')==0 && iF == 3 % check fMRI folder for ftap only
                    disp([use_dir ' is not an existing directory']) % message
                    use_dir = fullfile(subj_dir,fMRI_dirs{iF}); % generate directory in use
                    iScan = nScansExpected; % set iScan to max value within the iScan loop, so there are not duplicates of ftap PARRECs
                    
                    clear PAR_list REC_list
                    PAR_list = dir(fullfile(use_dir,'*.PAR')); PAR_list = AK_sortStruct(PAR_list,1,'SENSE_'); % generate structure containing information about mat files
                    REC_list = dir(fullfile(use_dir,'*.REC')); REC_list = AK_sortStruct(REC_list,1,'SENSE_'); % generate structure containing information about mat files
                    if ~isempty(PAR_list) && ~isempty(REC_list)
                        for iPR = 1:length(REC_list); % cycle through PARRECs
                            if any(AK_whichPattern(REC_list(iPR).name,dynStr));
                                 
                                NDARdata{NDARrow,1} = subj_dirs{iS}; % store subject
                                NDARdata{NDARrow,2} = fMRI_dirs{iF}; % store fMRI
                                NDARdata{NDARrow,3} = scan_dirs{iScan}; % store scan
                                NDARdata{NDARrow,4} = PAR_list(iPR).date; % store date
                                NDARdata{NDARrow,5} = PAR_list(iPR).name; % store PAR filename
                                NDARdata{NDARrow,6} = REC_list(iPR).name; % store REC filename
                                nScans(iS) = nScans(iS)+1; % document number of scans
                                NDARrow = NDARrow+1; % move to next row of NDARdata
                                
                                clear source check
                                source = fullfile(use_dir,PAR_list(iPR).name); % generate source filename
                                check = copyfile(source,destination_dir,'f'); % copy file from source to destination
                                if check == 0 % did file fail to copy
                                    disp(['**failed to copy file' source]) % message
                                end

                                clear source check
                                source = fullfile(use_dir,REC_list(iPR).name); % generate source filename
                                check = copyfile(source,destination_dir,'f'); % copy file from source to destination
                                if check == 0 % did file fail to copy
                                    disp(['**failed to copy file' source]) % message
                                end

                            end
                        end
                    else
                        disp(['either PAR or REC file is missing for ' subj_dirs{iS} fMRI_dirs{iF}]) % message
                    end
                elseif exist(use_dir,'dir')==0 
                    disp([use_dir ' is not an existing directory']) % message
                else
                    clear PAR_list REC_list
                    PAR_list = dir(fullfile(use_dir,'*.PAR')); PAR_list = AK_sortStruct(PAR_list,1,'SENSE_'); % generate structure containing information about mat files
                    REC_list = dir(fullfile(use_dir,'*.REC')); REC_list = AK_sortStruct(REC_list,1,'SENSE_'); % generate structure containing information about mat files
                    if ~isempty(PAR_list) && ~isempty(REC_list)
                        for iPR = 1:length(REC_list); % cycle through PARRECs
                            if any(AK_whichPattern(REC_list(iPR).name,dynStr));
                              
                                NDARdata{NDARrow,1} = subj_dirs{iS}; % store subject
                                NDARdata{NDARrow,2} = fMRI_dirs{iF}; % store fMRI
                                NDARdata{NDARrow,3} = scan_dirs{iScan}; % store scan
                                NDARdata{NDARrow,4} = PAR_list(iPR).date; % store date
                                NDARdata{NDARrow,5} = PAR_list(iPR).name; % store PAR filename
                                NDARdata{NDARrow,6} = REC_list(iPR).name; % store REC filename
                                nScans(iS) = nScans(iS)+1; % document number of scans
                                NDARrow = NDARrow+1; % move to next row of NDARdata
                                
                                clear source check
                                source = fullfile(use_dir,PAR_list(iPR).name); % generate source filename
                                check = copyfile(source,destination_dir,'f'); % copy file from source to destination
                                if check == 0 % did file fail to copy
                                    disp(['**failed to copy file' source]) % message
                                end

                                clear source check
                                source = fullfile(use_dir,REC_list(iPR).name); % generate source filename
                                check = copyfile(source,destination_dir,'f'); % copy file from source to destination
                                if check == 0 % did file fail to copy
                                    disp(['**failed to copy file' source]) % message
                                end
                                
                            end
                        end
                    else
                        disp(['either PAR or REC file is missing for ' subj_dirs{iS} fMRI_dirs{iF}]) % message
                    end
                end
                iScan = iScan+1; % advance counter
            end
        end
        
        for iM = 1:length(MRS_dirs); % cycle through MRS folders
            MRS_dir = fullfile(subj_dir,MRS_dirs{iM}); % generate directory
            if exist(MRS_dir,'dir')==0
                disp(['could not find MRS folder for' subj_dirs{iS} MRS_dirs{iM}]) % message
            elseif exist(MRS_dir,'dir')==7
                for iB = 1:length(subMRS_dirs); % cycle through subMRSs
                    use_dir = fullfile(MRS_dir,subMRS_dirs{iB}); % generate directory in use
                    if exist(use_dir,'dir')==0
                        disp(['could not find subMRS folder for' subj_dirs{iS} MRS_dirs{iM} subMRS_dirs{iB}]) % message
                    elseif exist(use_dir,'dir')==7 
                        clear actSPAR_list actSDAT_list refSPAR_list refSDAT_list
                        actSPAR_list = dir(fullfile(use_dir,'*act.SPAR')); % generate structure containing information about files
                        actSDAT_list = dir(fullfile(use_dir,'*act.SDAT')); % generate structure containing information about files
                        refSPAR_list = dir(fullfile(use_dir,'*ref.SPAR')); % generate structure containing information about files
                        refSDAT_list = dir(fullfile(use_dir,'*ref.SDAT')); % generate structure containing information about files
                        if ~isempty(actSPAR_list) && ~isempty(actSDAT_list)

                            NDARdata{NDARrow,1} = subj_dirs{iS}; % store subject
                            NDARdata{NDARrow,2} = MRS_dirs{iM}; % store MRS
                            NDARdata{NDARrow,3} = subMRS_dirs{iB}; % store subMRS
                            NDARdata{NDARrow,4} = actSPAR_list.date;
                            NDARdata{NDARrow,7} = actSPAR_list.name; % store actSPAR filename
                            NDARdata{NDARrow,8} = actSDAT_list.name; % store actSDAT filename
                            NDARdata{NDARrow,9} = refSPAR_list.name; % store refSPAR filename
                            NDARdata{NDARrow,10} = refSDAT_list.name; % store refSDAT filename
                            nMRSs(iS) = nMRSs(iS)+length(actSPAR_list); % document number of subMRSs
                            NDARrow = NDARrow+1; % move to next row of NDARdata

                            for iAP = 1:length(actSPAR_list); % cycle through PAR files
                                clear source check
                                source = fullfile(use_dir,actSPAR_list(iAP).name); % generate source filename
                                check = copyfile(source,destination_dir,'f'); % copy file from source to destination
                                if check == 0 % did file fail to copy
                                    disp(['**failed to copy file: ' source]) % message
                                end
                            end

                            for iAD = 1:length(actSDAT_list); % cycle through REC files
                                clear source check
                                source = fullfile(use_dir,actSDAT_list(iAD).name); % generate source filename
                                check = copyfile(source,destination_dir,'f'); % copy file from source to destination
                                if check == 0 % did file fail to copy
                                    disp(['**failed to copy file: ' source]) % message
                                end
                            end

                            for iRP = 1:length(refSPAR_list); % cycle through PAR files
                                clear source check
                                source = fullfile(use_dir,refSPAR_list(iRP).name); % generate source filename
                                check = copyfile(source,destination_dir,'f'); % copy file from source to destination
                                if check == 0 % did file fail to copy
                                    disp(['**failed to copy file: ' source]) % message
                                end
                            end

                            for iRD = 1:length(refSDAT_list); % cycle through REC files
                                clear source check
                                source = fullfile(use_dir,refSDAT_list(iRD).name); % generate source filename
                                check = copyfile(source,destination_dir,'f'); % copy file from source to destination
                                if check == 0 % did file fail to copy
                                    disp(['**failed to copy file: ' source]) % message
                                end
                            end

                        else
                            disp(['either SPAR or SDAT file is missing for ' subj_dirs{iS} MRS_dirs{iM}]) % message
                        end
                    end
                end
            end
        end
        
    end
end


end

