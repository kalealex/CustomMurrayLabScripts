% Create_fMRI_HeadMotion_Data_Table_GABAinASD
% Murray_Lab_2016
% Created by AMK on 7/12/16
% Modified 2/2/17

%% designate directory, subject, and logfile indentifying information

top_dir = 'L:\MurrayLab\ASD\Data';
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab';
% list current subjects here; won't work for subjects with unanalyzed fMRI data
% all subjects
subjects = dir(fullfile(top_dir,'G*')); % look for folders matching subject code format
subjects = {subjects.name}; % cat names
subjects = [subjects {'KG122'}]; % add kiddo
fmri_dirs = {'fMRI1','fMRI2','ftap'};

% strings for logfile/PRT file recognition and association
MTlocStr = {'MTlocalizer'};
V1locStr = {'V1localizer'};
V1_fixlocStr = {'V1localizer_fix'};
contrastStr = {'contrast'};
supStr = {'suppression'};
sumStr = {'summation'};
ftapStr = {'ftap'};

%% designate table structure and field names

% name columns of eye tracking data tables
table_columns = {'subject','session','fMRI set#','framewise displacement RMS','framewise displacement SD','max framewise displacement'};

% create cell arrays to store data and save into spreadsheets
MTloc = cell(3*length(subjects)+50,length(table_columns)); % 50 rows for buffer
V1loc = cell(3*length(subjects)+50,length(table_columns)); % 50 rows for buffer
V1_fixloc = cell(3*length(subjects)+50,length(table_columns)); % 50 rows for buffer
contrast = cell(6*length(subjects)+100,length(table_columns)); % 100 rows for buffer
sup = cell(6*length(subjects)+100,length(table_columns)); % 100 rows for buffer; sup and sum have the same condition names but separate cell arrays
sum = cell(6*length(subjects)+100,length(table_columns)); % 100 rows for buffer
ftap = cell(6*length(subjects)+150,length(table_columns)); % 150 rows for buffer

% add column names to cell arrays
MTloc(1,:) = table_columns;
V1loc(1,:) = table_columns;
V1_fixloc(1,:) = table_columns;
contrast(1,:) = table_columns;
sup(1,:) = table_columns;
sum(1,:) = table_columns;
ftap(1,:) = table_columns;

% create counters for rows of each cell array
MTlocRow = 2;
V1locRow = 2;
V1_fixlocRow = 2;
contrastRow = 2;
supRow = 2;
sumRow = 2;
ftapRow = 2;

%% get list of assicated prt filenames for each set

prtList = AK_GABAinASD_getAllAssociatedPRTnamesBySetSessionSubject(subjects,false);

%% pull data from data structure into appropriat spreadsheets

for iS = 1:length(subjects) % cycle through subjects 
    for iF = 1:length(fmri_dirs) % cycle through fMRI sessions
        if exist(fullfile(top_dir,subjects{iS},fmri_dirs{iF}),'file')
            % get head motion framewise displacement
            clear FD FDrms FDsd FDmax setlist emptyIndex
            [FD,FDrms,FDsd,FDmax,setlist] = AK_GABAinASD_getFDinfMRI(subjects{iS},fmri_dirs{iF});
            % remove nans and empty cells
            emptyIndex = cellfun(@isempty,FD);
            FD(emptyIndex) = [];
            FDrms(emptyIndex) = [];
            FDsd(emptyIndex) = [];
            FDmax(emptyIndex) = [];

            % create prtList for the current subject and session
            clear SubjPRTlist
            SubjPRTlist = prtList((AK_findStrMatch(prtList(:,1),subjects{iS}) & AK_findStrMatch(prtList(:,2),fmri_dirs{iF})),:);

            for iSet = 1:length(setlist) % cycle through asc files
                % find and store setN as a variable
                clear setNidx
                setNidx = strfind(setlist{iSet},'set');
                if ~isnan(setlist{iSet}(setNidx+length('set')))
                    setN = str2double(setlist{iSet}(setNidx+length('set')));
                end

                % find the SubjPRTlist row that corresponds to set iSet
                clear SubjPRTlistRow
                SubjPRTlistRow = cellfun(@(x) x==setN,SubjPRTlist(:,3),'UniformOutput',1);

                % determine which type of set iSet corresponds to
                if AK_whichPattern(SubjPRTlist{SubjPRTlistRow,4},MTlocStr,true)==1 % is this an MT localizer?

                    % checks to advance row counter and avoid overwriting
                    if all(~cellfun(@isempty,MTloc(MTlocRow,:))) % check to advance counter
                        MTlocRow = MTlocRow+1; % advance counter
                    elseif any(~cellfun(@isempty,MTloc(MTlocRow,:))) % check that data is not being overwritten
                        disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,4} ' at row ' num2str(MTlocRow) ' of MTloc']) % message
                        MTlocRow = MTlocRow+1; % advance counter
                    end

                    % store stats about eye tracking data
                    MTloc{MTlocRow,1} = subjects{iS}; % store subject name
                    MTloc{MTlocRow,2} = fmri_dirs{iF}; % store fMRI session name
                    MTloc{MTlocRow,3} = num2str(setN); % store set#
                    MTloc{MTlocRow,4} = FDrms(iSet);
                    MTloc{MTlocRow,5} = FDsd(iSet);
                    MTloc{MTlocRow,6} = FDmax(iSet);

                    % remove corresponding row from SubjPRTlist
                    SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

                elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,4},V1_fixlocStr,true)==1 % is this an V1_fix localizer?; must check before V1 localizer because of overlap in strings

                    % checks to advance row counter and avoid overwriting
                    if all(~cellfun(@isempty,V1_fixloc(V1_fixlocRow,:))) % check to advance counter
                        V1_fixlocRow = V1_fixlocRow+1; % advance counter
                    elseif any(~cellfun(@isempty,V1_fixloc(V1_fixlocRow,:))) % check that data is not being overwritten
                        disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,4} ' at row ' num2str(V1_fixlocRow) ' of V1_fixloc']) % message
                        V1_fixlocRow = V1_fixlocRow+1; % advance counter
                    end

                    % store stats about eye tracking data
                    V1_fixloc{V1_fixlocRow,1} = subjects{iS}; % store subject name
                    V1_fixloc{V1_fixlocRow,2} = fmri_dirs{iF}; % store fMRI session name
                    V1_fixloc{V1_fixlocRow,3} = num2str(setN); % store set#
                    V1_fixloc{V1_fixlocRow,4} = FDrms(iSet);
                    V1_fixloc{V1_fixlocRow,5} = FDsd(iSet);
                    V1_fixloc{V1_fixlocRow,6} = FDmax(iSet);

                    % remove corresponding row from SubjPRTlist
                    SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

                elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,4},V1locStr,true)==1 % is this an V1 localizer?

                    % checks to advance row counter and avoid overwriting
                    if all(~cellfun(@isempty,V1loc(V1locRow,:))) % check to advance counter
                        V1locRow = V1locRow+1; % advance counter
                    elseif any(~cellfun(@isempty,V1loc(V1locRow,:))) % check that data is not being overwritten
                        disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,4} ' at row ' num2str(V1locRow) ' of V1loc']) % message
                        V1locRow = V1locRow+1; % advance counter
                    end

                    % store stats about eye tracking data
                    V1loc{V1locRow,1} = subjects{iS}; % store subject name
                    V1loc{V1locRow,2} = fmri_dirs{iF}; % store fMRI session name
                    V1loc{V1locRow,3} = num2str(setN); % store set#
                    V1loc{V1locRow,4} = FDrms(iSet);
                    V1loc{V1locRow,5} = FDsd(iSet);
                    V1loc{V1locRow,6} = FDmax(iSet);

                    % remove corresponding row from SubjPRTlist
                    SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

                elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,4},contrastStr,true)==1 % is this a contrast scan?

                    % checks to advance row counter and avoid overwriting
                    if all(~cellfun(@isempty,contrast(contrastRow,:))) % check to advance counter
                        contrastRow = contrastRow+1; % advance counter
                    elseif any(~cellfun(@isempty,contrast(contrastRow,:))) % check that data is not being overwritten
                        disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,4} ' at row ' num2str(contrastRow) ' of contrast']) % message
                        contrastRow = contrastRow+1; % advance counter
                    end

                   % store stats about eye tracking data
                    contrast{contrastRow,1} = subjects{iS}; % store subject name
                    contrast{contrastRow,2} = fmri_dirs{iF}; % store fMRI session name
                    contrast{contrastRow,3} = num2str(setN); % store set#
                    contrast{contrastRow,4} = FDrms(iSet);
                    contrast{contrastRow,5} = FDsd(iSet);
                    contrast{contrastRow,6} = FDmax(iSet);

                    % remove corresponding row from SubjPRTlist
                    SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

                elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,4},supStr,true)==1 % is this a suppression scan?

                    % checks to advance row counter and avoid overwriting
                    if all(~cellfun(@isempty,sup(supRow,:))) % check to advance counter
                        supRow = supRow+1; % advance counter
                    elseif any(~cellfun(@isempty,sup(supRow,:))) % check that data is not being overwritten
                        disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,4} ' at row ' num2str(supRow) ' of sup']) % message
                        supRow = supRow+1; % advance counter
                    end

                    % store stats about eye tracking data
                    sup{supRow,1} = subjects{iS}; % store subject name
                    sup{supRow,2} = fmri_dirs{iF}; % store fMRI session name
                    sup{supRow,3} = num2str(setN); % store set#
                    sup{supRow,4} = FDrms(iSet);
                    sup{supRow,5} = FDsd(iSet);
                    sup{supRow,6} = FDmax(iSet);

                    % remove corresponding row from SubjPRTlist
                    SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

                elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,4},sumStr,true)==1 % is this a summation scan?

                    % checks to advance row counter and avoid overwriting
                    if all(~cellfun(@isempty,sum(sumRow,:))) % check to advance counter
                        sumRow = sumRow+1; % advance counter
                    elseif any(~cellfun(@isempty,sum(sumRow,:))) % check that data is not being overwritten
                        disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,4} ' at row ' num2str(sumRow) ' of sum']) % message
                        sumRow = sumRow+1; % advance counter
                    end

                    % store stats about eye tracking data
                    sum{sumRow,1} = subjects{iS}; % store subject name
                    sum{sumRow,2} = fmri_dirs{iF}; % store fMRI session name
                    sum{sumRow,3} = num2str(setN); % store set#
                    sum{sumRow,4} = FDrms(iSet);
                    sum{sumRow,5} = FDsd(iSet);
                    sum{sumRow,6} = FDmax(iSet);

                    % remove corresponding row from SubjPRTlist
                    SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist
                    
                elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,4},ftapStr,true)==1 % is this a ftap scan?

                    % checks to advance row counter and avoid overwriting
                    if all(~cellfun(@isempty,ftap(ftapRow,:))) % check to advance counter
                        ftapRow = ftapRow+1; % advance counter
                    elseif any(~cellfun(@isempty,ftap(ftapRow,:))) % check that data is not being overwritten
                        disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,4} ' at row ' num2str(ftapRow) ' of ftap']) % message
                        ftapRow = ftapRow+1; % advance counter
                    end

                    % store stats about eye tracking data
                    ftap{ftapRow,1} = subjects{iS}; % store subject name
                    ftap{ftapRow,2} = fmri_dirs{iF}; % store fMRI session name
                    ftap{ftapRow,3} = num2str(setN); % store set#
                    ftap{ftapRow,4} = FDrms(iSet);
                    ftap{ftapRow,5} = FDsd(iSet);
                    ftap{ftapRow,6} = FDmax(iSet);

                    % remove corresponding row from SubjPRTlist
                    SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist
                    
                end    
            end
            if ~isempty(SubjPRTlist)
                disp('***Some unmatched set#s remaining in prtList:'); % message
                SubjPRTlist
            end
        end
    end
end

% trim empty rows from end of each cell array
MTloc(all(cellfun(@isempty,MTloc(:,1:2)),2),:) = [];
V1loc(all(cellfun(@isempty,V1loc(:,1:2)),2),:) = [];
V1_fixloc(all(cellfun(@isempty,V1_fixloc(:,1:2)),2),:) = [];
contrast(all(cellfun(@isempty,contrast(:,1:2)),2),:) = [];
sup(all(cellfun(@isempty,sup(:,1:2)),2),:) = [];
sum(all(cellfun(@isempty,sum(:,1:2)),2),:) = [];
ftap(all(cellfun(@isempty,ftap(:,1:2)),2),:) = [];

%% write cell arrays to xls file spreadsheets

success = zeros(6,1); % preallocate
% xlsFilename = fullfile(top_dir,'fMRI_EyeTracking_Data_Tables.xlsx'); % name file
xlsFilename = fullfile(top_dir,'fMRI_HeadMotion_Data_Tables.xlsx'); % name file for only subjects with good data

success(1) = xlswrite(xlsFilename,MTloc,'MT Loc');
success(2) = xlswrite(xlsFilename,V1loc,'V1 Loc');
success(3) = xlswrite(xlsFilename,V1_fixloc,'V1_fix Loc');
success(4) = xlswrite(xlsFilename,contrast,'Contrast');
success(5) = xlswrite(xlsFilename,sup,'Suppression');
success(6) = xlswrite(xlsFilename,sum,'Summation');
success(7) = xlswrite(xlsFilename,ftap,'ftap');
