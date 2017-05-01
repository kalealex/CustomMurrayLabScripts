function AK_fMRImetadata_forNDAR( metadata_filename, subjInfo_filename )
%AK_fMRImetadata_forNDAR prepares a .csv file containing appropriately 
% formatted meta-data for all necessary fMRI data files and saves this file 
% to the NDAR submission folder with the fullfile name specified in the 
% first function input (metadata_filename). It does this for the subjects 
% whose information is contained in the file designated in the second 
% function input (subjInfo_filename). See AK_prepareNDARsubmission.m 
% documentation for more detailed instructions on the NDAR submission 
% process and file formatting.
%
% HISTORY:
% used to be fMRI_Data_for_NDAR
% Murray Lab 2016
% Created by Alex Kale 1/8/16

%% establish directory by subjects

top_dir = 'L:\MurrayLab\ASD\Data'; % set up base directory
fMRI_dirs = {'fMRI1','fMRI2','ftap'}; % designate folders for sessions
scan_dirs = {'3danat','set1','set2','set3','set4','set5','set6','set7','set8','set9'}; % designate folders for scans
% for filling experiment id, scan type, and study fields
associatedPRTstr = {'MTlocalizer*','contrast*','suppression*','summation*','V1localizer*','V1localizer_fix*'};
expIDs = {419,416,417,418,420,421}; % same indexing as scanStr (see below)
expTRs = {125,125,125,125,80,80}; 

studyStr = {'MT localizer','Contrast','Suppression','Summation','V1 localizer','V1 localizer-fixation'};
dynStr = {'MPRAGE','125Dyn','80Dyn','170Dyn'};

% subjInfo_filename = fullfile(top_dir,'NDARnotes','NDAR_SubjectInfo.xlsx'); % create full filename for NDAR subject info spreadsheet
[~,~,subjInfo] = xlsread(subjInfo_filename); % load subject info
subj_dirs = subjInfo(2:end,1); % designate subject ids

dateFormat = 'mm/dd/yyyy'; % string for date formatting

%% create cell array for data output

NDARdata = cell(length(subj_dirs)*length(fMRI_dirs)*length(scan_dirs)+1,70);
NDARdata(1,1:2) = {'image','3'};
NDARdata(2,:) = {'subjectkey','src_subject_id','interview_date','interview_age','gender','comments_misc','image_file','image_thumbnail_file','image_description','experiment_id','scan_type','scan_object','image_file_format','data_file2','data_file2_type','image_modality','scanner_manufacturer_pd','scanner_type_pd','scanner_software_versions_pd','magnetic_field_strength','mri_repetition_time_pd','mri_echo_time_pd','flip_angle','acquisition_matrix','mri_field_of_view_pd','patient_position','photomet_interpret','receive_coil','transmit_coil','transformation_performed','transformation_type','image_history','image_num_dimensions','image_extent1','image_extent2','image_extent3','image_extent4','extent4_type','image_extent5','extent5_type','image_unit1','image_unit2','image_unit3','image_unit4','image_unit5','image_resolution1','image_resolution2','image_resolution3','image_resolution4','image_resolution5','image_slice_thickness','image_orientation','qc_outcome','qc_description','qc_fail_quest_reason','decay_correction','frame_end_times','frame_end_unit','frame_start_times','frame_start_unit','pet_isotope','pet_tracer','time_diff_inject_to_image','time_diff_units','pulse_seq','slice_acquisition','software_preproc','study','week','experiment_description'};

%% load data and organize in table, per file, per subject

NDARrow = 3; % designate variable for indexing NDARdata rows

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
                                NDARdata{NDARrow,1} = subjInfo{iS+1,2}; % store GUID
                                NDARdata{NDARrow,2} = subj_dirs{iS}; % store subject id
                                NDARdata{NDARrow,3} = datestr(REC_list(iPR).date,dateFormat); % store interview date
                                NDARdata{NDARrow,4} = AK_calculateInterviewAge_forNDAR(subjInfo{iS+1,3}, datestr(REC_list(iPR).date,dateFormat)); % store age
                                NDARdata{NDARrow,5} = subjInfo{iS+1,4}; % store gender
                                NDARdata{NDARrow,7} = REC_list(iPR).name; % store REC filename (image file)
                                NDARdata{NDARrow,9} = 'fMRI'; % store image description
                                NDARdata{NDARrow,12} = 'Live'; % store scan object
                                NDARdata{NDARrow,13} = 'PARREC'; % store image file format
                                NDARdata{NDARrow,14} = PAR_list(iPR).name; % store PAR filename (data file2)
                                NDARdata{NDARrow,15} = 'Associated PAR file'; % store data file2 type
                                NDARdata{NDARrow,16} = 'MRI'; % store image modality
                                %%% SCANNER SPECS ADDED APRIL, 2017
                                NDARdata{NDARrow,17} = 'Philips'; % store scanner manufacturer
                                NDARdata{NDARrow,18} = 'Achieva'; % store scanner type
                                if datetime(PAR_list(iPR).date) < datetime('9/15/2015') % software changed from 3.2.3 to 5.1.7 on 9/15/2015
                                    NDARdata{NDARrow,19} = '3.2.3'; % store scanner software version
                                else
                                    NDARdata{NDARrow,19} = '5.1.7'; % store scanner software version
                                end
                                NDARdata{NDARrow,20} = '3T'; % store magnetic field strength
                                NDARdata{NDARrow,26} = 'head first supine'; % store patient position
                                NDARdata{NDARrow,27} = 'monochrome 2'; % store description of photometric interpretation
                                %%% 
                                NDARdata{NDARrow,28} = '32 channel headcoil'; % store receive coil
                                NDARdata{NDARrow,29} = 'headcoil'; % store transmit coil
                                NDARdata{NDARrow,30} = 'No'; % store transformation performed
                                NDARdata{NDARrow,53} = 'pass'; % store qc outcome
                                NDARdata{NDARrow,66} = 1; % store slice acquisition
                                NDARdata{NDARrow,67} = 'BrainVoyager QX'; % store software preproc

                                if iScan == 1 % if MPRAGE
                                    % store no experiment id
                                    NDARdata{NDARrow,11} = 'MR structural (MPRAGE)'; % store scan type
                                    %%% SCAN SPECS ADDED APRIL, 2017
                                    NDARdata{NDARrow,21} = 0.0076; % store repetition time (seconds)
                                    NDARdata{NDARrow,22} = 0.0035; % store echo time (seconds)
                                    NDARdata{NDARrow,23} = '7 degrees'; % store flip angle
                                    NDARdata{NDARrow,24} = '256 * 256 * 176 voxels'; % store acquisition matrix
                                    NDARdata{NDARrow,25} = '256 * 256 * 176 millimeters'; % store field of view
                                    %%%
                                    %%% INFORMATION ON IMAGE DIMENSIONS ADDED JAN, 2017
                                    NDARdata{NDARrow,33} = 3; % store number of image dimensions
                                    NDARdata{NDARrow,34} = 256; % store extent of image dimension 1
                                    NDARdata{NDARrow,35} = 256; % store extent of image dimension 2
                                    NDARdata{NDARrow,36} = 176; % store extent of image dimension 3
                                    % store no extent of image dimension 4
                                    % store no type for extent of image dimension 4 ('Time')
                                    NDARdata{NDARrow,41} = 'Millimeters'; % store units of image dimension 1
                                    NDARdata{NDARrow,42} = 'Millimeters'; % store units of image dimension 2
                                    NDARdata{NDARrow,43} = 'Millimeters'; % store units of image dimension 3
                                    % store no units of image dimension 4 ('Seconds')
                                    NDARdata{NDARrow,46} = 1; % store resolution of image dimension 1
                                    NDARdata{NDARrow,47} = 1; % store resolution of image dimension 2
                                    NDARdata{NDARrow,48} = 1; % store resolution of image dimension 3
                                    % store no resolution of image dimension 4
                                    NDARdata{NDARrow,51} = 1; % store image slice thickness
                                    NDARdata{NDARrow,52} = 'Sagittal'; % store image orientation
                                    %%%
                                    % store no study
                                else
                                    NDARdata{NDARrow,10} = 422; % store experiment id
                                    NDARdata{NDARrow,11} = 'MR structural (T2)'; % store scan type
                                    %%% SCAN SPECS ADDED APRIL, 2017
                                    NDARdata{NDARrow,21} = 2; % store repetition time (seconds)
                                    NDARdata{NDARrow,22} = 0.025; % store echo time (seconds)
                                    NDARdata{NDARrow,23} = '79 degrees'; % store flip angle
                                    NDARdata{NDARrow,24} = '80 * 80 * 30 voxels'; % store acquisition matrix
                                    NDARdata{NDARrow,25} = '240 * 240 * 105 millimeters'; % store field of view
                                    %%%
                                    %%% INFORMATION ON IMAGE DIMENSIONS ADDED JAN, 2017
                                    NDARdata{NDARrow,33} = 4; % store number of image dimensions
                                    NDARdata{NDARrow,34} = 80; % store extent of image dimension 1
                                    NDARdata{NDARrow,35} = 80; % store extent of image dimension 2
                                    NDARdata{NDARrow,36} = 30; % store extent of image dimension 3
                                    NDARdata{NDARrow,37} = 170; % store extent of image dimension 4
                                    NDARdata{NDARrow,38} = 'Time';% store type for extent of image dimension 4
                                    NDARdata{NDARrow,41} = 'Millimeters'; % store units of image dimension 1
                                    NDARdata{NDARrow,42} = 'Millimeters'; % store units of image dimension 2
                                    NDARdata{NDARrow,43} = 'Millimeters'; % store units of image dimension 3
                                    NDARdata{NDARrow,44} = 'Seconds';% store units of image dimension 4
                                    NDARdata{NDARrow,46} = 3; % store resolution of image dimension 1
                                    NDARdata{NDARrow,47} = 3; % store resolution of image dimension 2
                                    NDARdata{NDARrow,48} = 3.5; % store resolution of image dimension 3
                                    NDARdata{NDARrow,49} = 2;% store resolution of image dimension 4
                                    NDARdata{NDARrow,51} = 3; % store image slice thickness
                                    NDARdata{NDARrow,52} = 'Axial'; % store image orientation
                                    %%%
                                    NDARdata{NDARrow,68} = 'ftap'; % store study
                                end

                                NDARrow = NDARrow+1; % move to next row of NDARdata
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
                                NDARdata{NDARrow,1} = subjInfo{iS+1,2}; % store GUID
                                NDARdata{NDARrow,2} = subj_dirs{iS}; % store subject id
                                NDARdata{NDARrow,3} = datestr(REC_list(iPR).date,dateFormat); % store interview date
                                NDARdata{NDARrow,4} = AK_calculateInterviewAge_forNDAR(subjInfo{iS+1,3}, datestr(REC_list(iPR).date,dateFormat)); % store age
                                NDARdata{NDARrow,5} = subjInfo{iS+1,4}; % store gender
                                NDARdata{NDARrow,7} = REC_list(iPR).name; % store REC filename (image file)
                                NDARdata{NDARrow,9} = 'fMRI'; % store image description
                                NDARdata{NDARrow,12} = 'Live'; % store scan object
                                NDARdata{NDARrow,13} = 'PARREC'; % store image file format
                                NDARdata{NDARrow,14} = PAR_list(iPR).name; % store PAR filename (data file2)
                                NDARdata{NDARrow,15} = 'Associated PAR file'; % store data file2 type
                                NDARdata{NDARrow,16} = 'MRI'; % store image modality
                                %%% SCANNER SPECS ADDED APRIL, 2017
                                NDARdata{NDARrow,17} = 'Philips'; % store scanner manufacturer
                                NDARdata{NDARrow,18} = 'Achieva'; % store scanner type
                                if datetime(PAR_list(iPR).date) < datetime('9/15/2015') % software changed from 3.2.3 to 5.1.7 on 9/15/2015
                                    NDARdata{NDARrow,19} = '3.2.3'; % store scanner software version
                                else
                                    NDARdata{NDARrow,19} = '5.1.7'; % store scanner software version
                                end
                                NDARdata{NDARrow,20} = '3T'; % store magnetic field strength
                                NDARdata{NDARrow,26} = 'head first supine'; % store patient position
                                NDARdata{NDARrow,27} = 'monochrome 2'; % store description of photometric interpretation
                                %%% 
                                NDARdata{NDARrow,28} = '32 channel headcoil'; % store receive coil
                                NDARdata{NDARrow,29} = 'headcoil'; % store transmit coil
                                NDARdata{NDARrow,30} = 'No'; % store transformation performed
                                NDARdata{NDARrow,53} = 'pass'; % store qc outcome
                                NDARdata{NDARrow,66} = 1; % store slice acquisition
                                NDARdata{NDARrow,67} = 'BrainVoyager QX'; % store software preproc

                                if iF == 3; % if ftap
                                    if iScan == 1 % if MPRAGE
                                        % store no experiment id
                                        NDARdata{NDARrow,11} = 'MR structural (MPRAGE)'; % store scan type
                                        %%% SCAN SPECS ADDED APRIL, 2017
                                        NDARdata{NDARrow,21} = 0.0076; % store repetition time (seconds)
                                        NDARdata{NDARrow,22} = 0.0035; % store echo time (seconds)
                                        NDARdata{NDARrow,23} = '7 degrees'; % store flip angle
                                        NDARdata{NDARrow,24} = '256 * 256 * 176 voxels'; % store acquisition matrix
                                        NDARdata{NDARrow,25} = '256 * 256 * 176 millimeters'; % store field of view
                                        %%%
                                        %%% INFORMATION ON IMAGE DIMENSIONS ADDED JAN, 2017
                                        NDARdata{NDARrow,33} = 3; % store number of image dimensions
                                        NDARdata{NDARrow,34} = 256; % store extent of image dimension 1
                                        NDARdata{NDARrow,35} = 256; % store extent of image dimension 2
                                        NDARdata{NDARrow,36} = 176; % store extent of image dimension 3
                                        % store no extent of image dimension 4
                                        % store no type for extent of image dimension 4 ('Time')
                                        NDARdata{NDARrow,41} = 'Millimeters'; % store units of image dimension 1
                                        NDARdata{NDARrow,42} = 'Millimeters'; % store units of image dimension 2
                                        NDARdata{NDARrow,43} = 'Millimeters'; % store units of image dimension 3
                                        % store no units of image dimension 4 ('Seconds')
                                        NDARdata{NDARrow,46} = 1; % store resolution of image dimension 1
                                        NDARdata{NDARrow,47} = 1; % store resolution of image dimension 2
                                        NDARdata{NDARrow,48} = 1; % store resolution of image dimension 3
                                        % store no resolution of image dimension 4
                                        NDARdata{NDARrow,51} = 1; % store image slice thickness
                                        NDARdata{NDARrow,52} = 'Sagittal'; % store image orientation
                                        %%%
                                        % store no study
                                    else
                                        NDARdata{NDARrow,10} = 422; % store experiment id
                                        NDARdata{NDARrow,11} = 'MR structural (T2)'; % store scan type
                                        %%% SCAN SPECS ADDED APRIL, 2017
                                        NDARdata{NDARrow,21} = 2; % store repetition time (seconds)
                                        NDARdata{NDARrow,22} = 0.025; % store echo time (seconds)
                                        NDARdata{NDARrow,23} = '79 degrees'; % store flip angle
                                        NDARdata{NDARrow,24} = '80 * 80 * 30 voxels'; % store acquisition matrix
                                        NDARdata{NDARrow,25} = '240 * 240 * 105 millimeters'; % store field of view
                                        %%%
                                        %%% INFORMATION ON IMAGE DIMENSIONS ADDED JAN, 2017
                                        NDARdata{NDARrow,33} = 4; % store number of image dimensions
                                        NDARdata{NDARrow,34} = 80; % store extent of image dimension 1
                                        NDARdata{NDARrow,35} = 80; % store extent of image dimension 2
                                        NDARdata{NDARrow,36} = 30; % store extent of image dimension 3
                                        NDARdata{NDARrow,37} = 170; % store extent of image dimension 4
                                        NDARdata{NDARrow,38} = 'Time';% store type for extent of image dimension 4
                                        NDARdata{NDARrow,41} = 'Millimeters'; % store units of image dimension 1
                                        NDARdata{NDARrow,42} = 'Millimeters'; % store units of image dimension 2
                                        NDARdata{NDARrow,43} = 'Millimeters'; % store units of image dimension 3
                                        NDARdata{NDARrow,44} = 'Seconds';% store units of image dimension 4
                                        NDARdata{NDARrow,46} = 3; % store resolution of image dimension 1
                                        NDARdata{NDARrow,47} = 3; % store resolution of image dimension 2
                                        NDARdata{NDARrow,48} = 3.5; % store resolution of image dimension 3
                                        NDARdata{NDARrow,49} = 2;% store resolution of image dimension 4
                                        NDARdata{NDARrow,51} = 3; % store image slice thickness
                                        NDARdata{NDARrow,52} = 'Axial'; % store image orientation
                                        %%%
                                        NDARdata{NDARrow,68} = 'ftap'; % store study
                                    end
                                else % if fMRI 1 or 2
                                    if iScan == 1 % if MPRAGE
                                        % store no experiment id
                                        NDARdata{NDARrow,11} = 'MR structural (MPRAGE)'; % store scan type
                                        %%% SCAN SPECS ADDED APRIL, 2017
                                        NDARdata{NDARrow,21} = 0.0076; % store repetition time (seconds)
                                        NDARdata{NDARrow,22} = 0.0035; % store echo time (seconds)
                                        NDARdata{NDARrow,23} = '7 degrees'; % store flip angle
                                        NDARdata{NDARrow,24} = '256 * 256 * 176 voxels'; % store acquisition matrix
                                        NDARdata{NDARrow,25} = '256 * 256 * 176 millimeters'; % store field of view
                                        %%%
                                        %%% INFORMATION ON IMAGE DIMENSIONS ADDED JAN, 2017
                                        NDARdata{NDARrow,33} = 3; % store number of image dimensions
                                        NDARdata{NDARrow,34} = 256; % store extent of image dimension 1
                                        NDARdata{NDARrow,35} = 256; % store extent of image dimension 2
                                        NDARdata{NDARrow,36} = 176; % store extent of image dimension 3
                                        % store no extent of image dimension 4
                                        % store no type for extent of image dimension 4 ('Time')
                                        NDARdata{NDARrow,41} = 'Millimeters'; % store units of image dimension 1
                                        NDARdata{NDARrow,42} = 'Millimeters'; % store units of image dimension 2
                                        NDARdata{NDARrow,43} = 'Millimeters'; % store units of image dimension 3
                                        % store no units of image dimension 4 ('Seconds')
                                        NDARdata{NDARrow,46} = 1; % store resolution of image dimension 1
                                        NDARdata{NDARrow,47} = 1; % store resolution of image dimension 2
                                        NDARdata{NDARrow,48} = 1; % store resolution of image dimension 3
                                        % store no resolution of image dimension 4
                                        NDARdata{NDARrow,51} = 1; % store image slice thickness
                                        NDARdata{NDARrow,52} = 'Sagittal'; % store image orientation
                                        %%%
                                        % store no study
                                    else
                                        try
                                            clear setIDfilename PRTfilename
                                            setIDfilename = [fullfile(use_dir,scan_dirs{iScan}) '_scanID.mat'];
                                            if exist(setIDfilename,'file')==0; % does the set scanID file exist?
                                                clear vtcfileName vtc 
                                                BVQXfile(0, 'clearallobjects');
                                                vtcfileName = [fullfile(use_dir,scan_dirs{iScan}) '.vtc']; % name vtc files 
                    %                             vtc = BVQXfile(deblank(vtcfileName(files,:))); % fill out vtc structure for purpose of naming scans
                                                vtc = BVQXfile(deblank(vtcfileName)); % fill out vtc structure for purpose of naming scans
                                                PRTfilename = vtc.NameOfLinkedPRT;
                                                save(setIDfilename,'PRTfilename'); % save associated .prt filename
                                            else
                                                load(setIDfilename,'PRTfilename'); % load associated .prt filename
                                            end
                                        catch
                                            PRTfilename = 'nothing';
                                            disp(['trouble processing ' vtcfileName]); % message
                                        end
                                        for iType = 1:length(associatedPRTstr); % cycle through types of scan
                                            if ~isempty(regexp(PRTfilename,regexptranslate('wildcard',associatedPRTstr{iType}),'once')) % is scan of type associatedPRTstr{iType}?
                                                NDARdata{NDARrow,10} = expIDs{iType}; % store experiment id   DEPENDS ON TYPE OF SCAN
                                                %%% INFORMATION ON IMAGE DIMENSIONS ADDED JAN, 2017
                                                NDARdata{NDARrow,37} = expTRs{iType}; % store extent of image dimension 4   DEPENDS ON TYPE OF SCAN
                                                %%%
                                                NDARdata{NDARrow,68} = studyStr{iType}; % store study   DEPENDS ON TYPE OF SCAN
                                            end
                                        end
                                        NDARdata{NDARrow,11} = 'MR structural (T2)'; % store scan type
                                        %%% SCAN SPECS ADDED APRIL, 2017
                                        NDARdata{NDARrow,21} = 2; % store repetition time (seconds)
                                        NDARdata{NDARrow,22} = 0.025; % store echo time (seconds)
                                        NDARdata{NDARrow,23} = '79 degrees'; % store flip angle
                                        NDARdata{NDARrow,24} = '80 * 80 * 30 voxels'; % store acquisition matrix
                                        NDARdata{NDARrow,25} = '240 * 240 * 105 millimeters'; % store field of view
                                        %%%
                                        %%% INFORMATION ON IMAGE DIMENSIONS ADDED JAN, 2017
                                        NDARdata{NDARrow,33} = 4; % store number of image dimensions
                                        NDARdata{NDARrow,34} = 80; % store extent of image dimension 1
                                        NDARdata{NDARrow,35} = 80; % store extent of image dimension 2
                                        NDARdata{NDARrow,36} = 30; % store extent of image dimension 3
                                        NDARdata{NDARrow,38} = 'Time';% store type for extent of image dimension 4
                                        NDARdata{NDARrow,41} = 'Millimeters'; % store units of image dimension 1
                                        NDARdata{NDARrow,42} = 'Millimeters'; % store units of image dimension 2
                                        NDARdata{NDARrow,43} = 'Millimeters'; % store units of image dimension 3
                                        NDARdata{NDARrow,44} = 'Seconds';% store units of image dimension 4
                                        NDARdata{NDARrow,46} = 3; % store resolution of image dimension 1
                                        NDARdata{NDARrow,47} = 3; % store resolution of image dimension 2
                                        NDARdata{NDARrow,48} = 3.5; % store resolution of image dimension 3
                                        NDARdata{NDARrow,49} = 2;% store resolution of image dimension 4
                                        NDARdata{NDARrow,51} = 3; % store image slice thickness
                                        NDARdata{NDARrow,52} = 'Axial'; % store image orientation
                                        %%%
                                    end
                                end

                                NDARrow = NDARrow+1; % move to next row of NDARdata
                            end
                        end
                    else
                        disp(['either PAR or REC file is missing for ' subj_dirs{iS} fMRI_dirs{iF}]) % message
                    end
                end
                iScan = iScan+1; % advance counter
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

