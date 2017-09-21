% makeEyetrackingQualityTable_GABAinASD
% MurrayLab 2016
% created by AMK 9/7/16

% THIS CODE IS MEANT TO BE RUN ONE SECTION AT A TIME

% The TOP SECTION either creates a eyetracking data quality table or
% appends to the existing one, depending on the value of newList (true -->
% new; false --> append). Each set per fMRI session per subject gets a row
% in the eyetracking data quality table if the following conditions are
% met:
    % the fMRI data has been analyzed in Brain Voyager, resulting in subfolders named set1, set2, set3, etc.
    % the subjects are accounted for in the most recent version of the file
        % AllAssociatedPRTfilenames.mat (this is taken care of within this
        % script by the function call on
        % AK_GABAinASD_getAllAssociatedPRTnamesBySetSessionSubject.m); it
        % may be the case that empty cells in the cell array prtList need
        % to be filled in from the notes, otherwise the script will likely
        % throw an indexing error
    % the edf associated with the set have been converted into asc files using the EDF2ASC converter (from SR Research)
% The quality table produced by the first section of code is a template and
% needs to be filled in with data quality booleans based on notes taken
% during scans. Copying and pasting booleans from the previously saved
% version of the quality table can save a lot of time here.

% NOTE: Where there are multiple .log files for a single set, as is often
% the case when a Presentation scenario is restarted or rerun (i.e., G127
% ftap2), the sting "BAD" must be appended to the end of the .log file
% names. Otherwise, the script will throw an error.
    
% The SECOND SECTION contains a few lines of code which are meant to assist
% in manual entry of 0s and 1s in the last two columns of the eyetracking
% data quality table. The idea here is to extract the information saved in
% the notes .txt files in order to inform which eyetracking data is worth
% analyzing (isData = 1) vs data is likely dominated by noise (isData = 0)
% as well as which eyetracking data is well calibrated (isGoodData = 1) vs
% data which is mediocre at best (isGoodData = 0). This rating system is
% reiterated in the key below.

% The THIRD SECTION is used to save the qualityTable as well as a list of
% disqualified subjects.

% Key to qualityTable:
% no row exists for sets without corresponding edfs
% every set with an existing edf is coded (in the finished table)
% isData = 0 if the data collected is pure noise
% isData = 1 if there is meaningful data, even when calibration is poor or mediocre
% isGoodData = 1 only if data is meaningful and well calibrated (according to notes)

%% create tables

addpath(genpath('L:\MurrayLab\ASD\Data'));

% create new table or append to existing one
newList = true;
        
% dir
top_dir = 'L:\MurrayLab\ASD\Data';
tableFilename = fullfile(top_dir,'EyetrackingQuality.mat');
oldTableFilename = fullfile(top_dir,['EyetrackingQuality_' datestr(now,'yymmdd') '.mat']);

% subjects, sessions, sets
subjects = dir(fullfile(top_dir,'G*')); % look for folders matching subject code format
subjects = {subjects.name}; % cat names
subjects = [subjects {'KG122'}]; % add kiddo
sessions = {'fMRI1','fMRI2','ftap'};
sets = {'set1','set2','set3','set4','set5','set6','set7','set8','set9'};

% scan identifying string for logfile/PRT file recognition and association
greencircIDstr = {'MTlocalizer','V1localizer_fix','V1localizer','contrast','suppression','summation'}; % for greencirc
ftapIDstr = {'fingertapping3_1','fingertapping3_2','fingertapping3_3';'ftap1','ftap2','ftap3'}; % for ftap

% load existing table and append, or create new:
if newList % prepare a new list
    % load old list and save with today's data in order to avoid overwriting old manual entry
    load(tableFilename,'dqSubjects','qualityTable'); % load existing table
    qt = qualityTable; % store the old quality table in order to allow user to copy and paste old manually entered values into new table
    save(oldTableFilename,'dqSubjects','qualityTable');
    % prepare cell array (qualityTable)
    qualityTable = cell(length(subjects)*length(sessions)*length(sets)+1,7);
    qualityTable(1,:) = {'subject','session','set','ascfile','logfile','isData','isGoodData'}; % column labels
    
    qRow = 2; % counter for quality table rows
else
    if exist(tableFilename,'file') % check that file exists
        load(tableFilename,'dqSubjects','qualityTable'); % load existing table
        
        % find list of subjects not already included in the cell array
        qualityTableSubjects = unique(prtList(2:end,1)); % list of unique subject names in existing list
        subjects = setdiff(subjects,qualityTableSubjects); % make subj_dirs = list of previously unincluded subjects

        qRow = length(qualityTable(:,1))+1; % create row counter  
    else
        disp('***Could not find prtList file; creating new list') % message
        % prepare cell array (qualityTable
        qualityTable = cell(length(subjects)*length(sessions)*length(sets)+1,7);
        qualityTable(1,:) = {'subject','session','set','ascfile','logfile','isData','isGoodData'}; % column labels

        qRow = 2; % counter for quality table rows
    end
end

% dqSubjects (might need to update this)
dqSubjects = {'G103','G117','G118','G134','G310','G326','G344'};

% list of confusing logfiles and which index is correct
confusingLogfilesTable = {'G101_4*.log', 2;...      % add to this table to avoid manual input each time the script gets confused about which version of a logfile to use
                        'G104_1*.log', 2;...        % often times the right .log file is the complete one (the one with larger file size)
                        'G106_1*.log', 1;...        % if there is any ambiguity, check the notes .txt file for the subject in question
                        'G107_1*.log', 1;...
                        'G112_1*.log', 2;...
                        'G112_19*.log', 2;...
                        'G123_1*.log', 3;...
                        'G123_2*.log', 2;...
                        'G123_8*.log', 3;...
                        'G307_10*.log', 2;...
                        'G312_1*.log', 2;...
                        'G316_1*.log', 2;...
                        'G322_1*.log', 3;...
                        'G332_10*.log', 2;...
                        'G340_15*.log', 2;...
                        'G345_1*.log', 1;...
                        'G348_10*.log', 2;...
                        'G348_2*.log', 2};
confusingLogfiles = confusingLogfilesTable(:,1);
confusingLogfilesIdx = [confusingLogfilesTable{:,2}];

% get associated prt filename list for each set for these subjects
prtList = AK_GABAinASD_getAllAssociatedPRTnamesBySetSessionSubject(subjects, 0); % second argument should be 1 unless you are sure that prtList is up to date

% fill in qualityTable
for iSubj = 1:length(subjects)
    % string to be used to find scan# and ascfile, logfile names
    clear scanNfindStr
    scanNfindStr = [subjects{iSubj} '_']; 
    
    for iSess = 1:length(sessions)
        % create prtList for this subject
        clear SubjPRTlist
        SubjPRTlist = prtList((AK_findStrMatch(prtList(:,1),subjects{iSubj}) & AK_findStrMatch(prtList(:,2),sessions{iSess})),:); 
        
        % preallocate
        clear setAscLogMatch
        setAscLogMatch = cell(0);
        salmRow =1; % reset counter
        % only bother trying to match up sets to .asc files to .log files
        % if the .prt file links hav been generated
        if ~isempty(SubjPRTlist) && ~all(cellfun(@isempty,SubjPRTlist(:,4)))
            % load and sort lists of ascfiles
            clear asc_list
            asc_list = dir(fullfile(top_dir,subjects{iSubj},sessions{iSess},'log_files','*.asc')); % create list of ascfiles
            asc_list = AK_sortStruct(asc_list,1,scanNfindStr); % sort ascfiles by scan number

            for iA = 1:length(asc_list) % cycle through ascfiles 
                % locate strings used to find scan# in asc filename and identify BAD eyetracking, respectively
                clear scanNindex BADindex
                scanNindex = strfind(asc_list(iA).name,scanNfindStr);
                BADindex = strfind(asc_list(iA).name,'BAD');
                if ~isempty(scanNindex) && isempty(BADindex) % filter asc files by name
                    %find and store scan# as variable scanN
                    clear scanN
                    if ~isnan(str2double(asc_list(iA).name(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1))) % check for two digit scan#
                        scanN = str2double(asc_list(iA).name(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1));
                    elseif ~isnan(str2double(asc_list(iA).name(scanNindex+length(scanNfindStr)))) % check for one digit scan#
                        scanN = str2double(asc_list(iA).name(scanNindex+length(scanNfindStr)));
                    end
                    % find associated logfile name in order to name conditions
                    clear LogfilefindStr log_list logfileName missingLogfiles missingLogfileIndex
                    LogfilefindStr = [scanNfindStr num2str(scanN) '*.log']; % sting to find associated logfile
                    log_list = dir(fullfile(top_dir,subjects{iSubj},sessions{iSess},'log_files',LogfilefindStr)); % create list of logfiles which include LogfilefindStr
                    log_list = AK_sortStruct(log_list,1,scanNfindStr); % sort logfile matches by scan number
                    if length(log_list)>1 % check number of logfiles
                        if any(strcmp(confusingLogfiles,LogfilefindStr)) % try to avoid asking for user input by using hard-coded table of logfiles which get confused by this script and the approriate indices
                            selectedLogfileIndex = confusingLogfilesIdx(strcmp(confusingLogfiles,LogfilefindStr));
                            logfileName = log_list(selectedLogfileIndex).name;
                        else % otherwise, ask user which logfile to use
                            disp(['Found multiple logfiles that match ' LogfilefindStr ':']) % message
                            for iL = 1:length(log_list) % cycle through logfiles
                                disp([num2str(iL) ': ' log_list(iL).name]) % list logfile names by index number
                            end
                            selectedLogfileIndex = input('Enter index (as a double) to select logfile to use from the above list:'); % prompt input for correct logfile
                            logfileName = log_list(selectedLogfileIndex).name;
                        end
                    elseif length(log_list)==1
                        logfileName = log_list(1).name;
                    elseif isempty(log_list) % check missing logfile list
                        disp(['***Checking missing logfile list for associated logfile for ' asc_list(iA).name]) % message
                        clear missingLogfiles missingLogfileIndex
                        load('L:\MurrayLab\ASD\Data\missing_logfile_list.mat','missingLogfiles');
                        missingLogfileIndex = AK_findStrMatch(missingLogfiles,LogfilefindStr); % is the associated logfile name on the list of missing logfiles?
                        logfileName = missingLogfiles{missingLogfileIndex}; % assign missing logfilename to variable in order to determine which table to add eye tracking data to
                    else
                        disp(['***Could not find associated logfile for ' asc_list(iA).name]) % message
                        logfileName = [];
                    end
                    % use logfile name associated with each ascfile to find the set number in the prtfile name list:
                    clear scanIDmatchIdx prtListScanIdx setN
                    if ~isempty(logfileName) && any(cellfun(@(x) any(strfind(logfileName,x)),greencircIDstr)) % for greencirc
                        
                        % which kind of scan is this logfile for?
                        scanIDmatchIdx = find(cellfun(@(x) any(strfind(logfileName,x)),greencircIDstr),1);
                        % what is the associated set#?
                        % THIS SECTION OF CODE ASSUMES THAT EACH TYPE OF
                        % SCAN WAS RUN A CORRECT NUMBER OF TIMES; IF THIS
                        % IS NOT THE CASE, IT IS LIKELY THAT THE fMRI
                        % ANALYSIS HAS INCORRECTLY ASSIGNED A .PRT FILE
                        % (RESULTING IN A BIG IN THIS CODE); TO AVOID THIS,
                        % WRITE THE MATCHING .PRT FILE NAME FOR THE .LOG
                        % FILE INTO 'prtList' FOR NOW AND RERUN THE fMRI
                        % ANALYSIS WITH THE CORRECT LINK BETWEEN THE .VTC
                        % AND THE .PRT FILES
                        prtListScanIdx = find(AK_findStrMatch(SubjPRTlist(:,4),greencircIDstr{scanIDmatchIdx},1)); % create index to find set#
                        setN = SubjPRTlist{prtListScanIdx,3}; % store set# as variable
                        SubjPRTlist(prtListScanIdx(1),:) = [];% remove row from SubjPRTlist
                    elseif ~isempty(logfileName) && any(cellfun(@(x) any(strfind(logfileName,x)),ftapIDstr(1,:))) % for ftap (logfile and prt file names are different) 
                        % which kind of scan is this logfile for?
                        scanIDmatchIdx = find(cellfun(@(x) any(strfind(logfileName,x)),ftapIDstr(1,:)),1);
                        % what is the associated set#?
                        % THIS SECTION OF CODE ASSUMES THAT EACH TYPE OF
                        % SCAN WAS RUN A CORRECT NUMBER OF TIMES; IF THIS
                        % IS NOT THE CASE, IT IS LIKELY THAT THE fMRI
                        % ANALYSIS HAS INCORRECTLY ASSIGNED A .PRT FILE
                        % (RESULTING IN A BIG IN THIS CODE); TO AVOID THIS,
                        % WRITE THE MATCHING .PRT FILE NAME FOR THE .LOG
                        % FILE INTO 'prtList' FOR NOW AND RERUN THE fMRI
                        % ANALYSIS WITH THE CORRECT LINK BETWEEN THE .VTC
                        % AND THE .PRT FILES
                        prtListScanIdx = find(AK_findStrMatch(SubjPRTlist(:,4),ftapIDstr{2,scanIDmatchIdx},1)); % create index to find set#
                        setN = SubjPRTlist{prtListScanIdx,3}; % store set# as variable
                        SubjPRTlist(prtListScanIdx(1),:) = [];% remove row from SubjPRTlist
                    else 
                        setN = [];
                    end
                    % store matching information gleaned
                    setAscLogMatch(salmRow,1) = {['set' num2str(setN)]}; % setN
                    setAscLogMatch(salmRow,2) = {asc_list(iA).name}; % ascfile name
                    setAscLogMatch(salmRow,3) = {logfileName}; % logfile name
                    salmRow = salmRow+1; % advance counter
                end
            end
        end
        
        for iSet = 1:length(sets)
            if ~isempty(setAscLogMatch) && ~any(cellfun(@isempty,setAscLogMatch(:,1))) && any(ismember(sets(iSet),setAscLogMatch(:,1))) % only add rows for existing sets; this assumes prtList is complete
                % fill in first 3 columns
                qualityTable(qRow,1) = subjects(iSubj);
                qualityTable(qRow,2) = sessions(iSess);
                qualityTable(qRow,3) = sets(iSet);
                % use set# to find ascfile and logfile names
                clear setMatchIdx selectedSALMidx
                if any(AK_findStrMatch(setAscLogMatch(:,1),sets{iSet}))
                    setMatchIdx = find(AK_findStrMatch(setAscLogMatch(:,1),sets{iSet}),1);
                    qualityTable(qRow,4) = setAscLogMatch(setMatchIdx,2); % store ascfile name
                    qualityTable(qRow,5) = setAscLogMatch(setMatchIdx,3); % store logfile name
                else 
                    disp(['*** Could not find set#, ascfile, logfile match for ' subjects{iSubj} ' ' sessions{iSess} ' ' sets{iSet}]); % message
                    for iSALM = 1:length(setAscLogMatch(:,1)) % cycle through logfiles
                        disp([num2str(iSALM) ': ' setAscLogMatch{iSALM,1} ' ' setAscLogMatch{iSALM,2} ' ' setAscLogMatch{iSALM,3}]) % list setAscLogMatch table by index number
                    end
                    selectedSALMidx = input('Enter index (as a double) to select a row to use from the above list:'); % prompt input for correct table row
                    qualityTable(qRow,4) = setAscLogMatch(selectedSALMidx,2); % store ascfile name
                    qualityTable(qRow,5) = setAscLogMatch(selectedSALMidx,3); % store logfile name
                end
                % fill in data quality ratings to the extent that it is
                % possible to automate this part:
                clear oldQtRowIdx
                % find the row in the old qualityTable (qt), if any, which
                % corresponds to the current row in the current qualityTable
                oldQtRowIdx = arrayfun(@(x) isequal(qualityTable(qRow,1:5),qt(x,1:5)), 1:length(qt(:,1)));
                if any(oldQtRowIdx) && sum(oldQtRowIdx) == 1 % there is one and only one match for this row
                    % use old quality table to fill in data quality rating
                    % columns in new quality table
                    qualityTable(qRow,6:7) = qt(oldQtRowIdx,6:7);
                end
                qRow = qRow+1; % advance counter 
            end
        end
    end
    
end

% message to user
disp('The last two columns must be filled in manually from notes.'); % message
disp('Values filled in by this section of code reflect data quality rating filled into prior versions of the quality table.')
disp('The remaining values left blank should be filled in from the notes.'); 
disp('Where set directories do not exist rows must be added to the table and filled in manually.')

%% fill in cells (to make manual entry faster)

fillRows = 939:947;
fillColumns = 6:7;
fillValue = 1;

qualityTable(fillRows,fillColumns) = {fillValue};

% no row exists for sets without corresponding edfs
% isData = 0 if the data collected is pure noise
% isData = 1 if there is meaningful data, even when calibration is poor or mediocre
% isGoodData = 1 only if data is meaningful and well calibrated (according to notes)

%% save

% remove empty rows
emptyRows = find(arrayfun(@(x) all(cellfun(@isempty, qualityTable(x,:))), 1:length(qualityTable(:,1))));
qualityTable(emptyRows,:) = [];

save(fullfile(top_dir,'EyetrackingQuality.mat'),'dqSubjects','qualityTable');
