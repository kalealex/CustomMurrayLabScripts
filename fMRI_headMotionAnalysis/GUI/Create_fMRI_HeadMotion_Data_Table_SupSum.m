% Create_fMRI_HeadMotion_Data_Table_SupSum
% Murray_Lab_2016
% Created by AMK on 8/1/16

%% designate directory, subject, and logfile indentifying information

top_dir = 'L:\MurrayLab\ASD\SuppressionSummation';
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab';
save_dir = 'L:\MurrayLab\DataTablesForGUI';
% list current subjects here; won't work for subjects with unanalyzed fMRI data
% all subjects
subjects = {'S_AMK_20160401','S_AVF_20160427','S_AVF_20160510','S_AW_20160420','S_BK_20160517','S_DP_20160427','S_JCM_20160426','S_MN_20160523','S_MPS_20160325','S_RM_20160414','S_SOM_20160414'};

% strings for logfile/PRT file recognition and association
MTlocStr = {'MTlocalizer'};
V1locStr = {'V1localizer'};
V1loc_meridianStr = {'V1localizer_meridian'};
supStr = {'suppression'};
sumStr = {'summation'};

%% designate table structure and field names

% name columns of eye tracking data tables
table_columns = {'subject','fMRI set#','framewise displacement RMS','framewise displacement SD','max framewise displacement'};

% create cell arrays to store data and save into spreadsheets
MTloc = cell(3*length(subjects)+50,length(table_columns)); % 50 rows for buffer
V1loc = cell(3*length(subjects)+50,length(table_columns)); % 50 rows for buffer
V1loc_meridian = cell(3*length(subjects)+50,length(table_columns)); % 50 rows for buffer
sup = cell(6*length(subjects)+100,length(table_columns)); % 100 rows for buffer; sup and sum have the same condition names but separate cell arrays
sum = cell(6*length(subjects)+100,length(table_columns)); % 100 rows for buffer

% add column names to cell arrays
MTloc(1,:) = table_columns;
V1loc(1,:) = table_columns;
V1loc_meridian(1,:) = table_columns;
sup(1,:) = table_columns;
sum(1,:) = table_columns;

% create counters for rows of each cell array
MTlocRow = 2;
V1locRow = 2;
V1loc_meridianRow = 2;
supRow = 2;
sumRow = 2;

%% get list of assicated prt filenames for each set

prtList = AK_SupSum_getAllAssociatedPRTnamesBySetSubject(subjects,0);

%% pull data from data structure into appropriat spreadsheets

for iS = 1:length(subjects) % cycle through subjects 
    % get head motion framewise displacement
    clear FD FDrms FDsd FDmax setlist emptyIndex
    [FD,FDrms,FDsd,FDmax,setlist] = AK_SupSum_getFDinfMRI(subjects{iS});
    % remove nans and empty cells
    emptyIndex = cellfun(@isempty,FD);
    FD(emptyIndex) = [];
    FDrms(emptyIndex) = [];
    FDsd(emptyIndex) = [];
    FDmax(emptyIndex) = [];

    % create prtList for the current subject
    clear SubjPRTlist
    SubjPRTlist = prtList(AK_findStrMatch(prtList(:,1),subjects{iS}),:);

    for iSet = 1:length(setlist) % cycle through asc files
        % find and store setN as a variable
        clear setNidx
        setNidx = strfind(setlist{iSet},'set');
        if ~isnan(setlist{iSet}(setNidx+length('set')))
            setN = str2double(setlist{iSet}(setNidx+length('set')));
        end

        % find the SubjPRTlist row that corresponds to set iSet
        clear SubjPRTlistRow
        SubjPRTlistRow = cellfun(@(x) x==setN,SubjPRTlist(:,2),'UniformOutput',1);

        % determine which type of set iSet corresponds to
        if AK_whichPattern(SubjPRTlist{SubjPRTlistRow,3},MTlocStr)==1 % is this an MT localizer?

            % checks to advance row counter and avoid overwriting
            if all(~cellfun(@isempty,MTloc(MTlocRow,:))) % check to advance counter
                MTlocRow = MTlocRow+1; % advance counter
            elseif any(~cellfun(@isempty,MTloc(MTlocRow,:))) % check that data is not being overwritten
                disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,3} ' at row ' num2str(MTlocRow) ' of MTloc']) % message
                MTlocRow = MTlocRow+1; % advance counter
            end

            % store stats about eye tracking data
            MTloc{MTlocRow,1} = subjects{iS}; % store subject name
            MTloc{MTlocRow,2} = num2str(setN); % store set#
            MTloc{MTlocRow,3} = FDrms(iSet);
            MTloc{MTlocRow,4} = FDsd(iSet);
            MTloc{MTlocRow,5} = FDmax(iSet);

            % remove corresponding row from SubjPRTlist
            SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

        elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,3},V1loc_meridianStr)==1 % is this a V1_meridian localizer?; must check before V1 localizer because of overlap in strings

            % checks to advance row counter and avoid overwriting
            if all(~cellfun(@isempty,V1loc_meridian(V1loc_meridianRow,:))) % check to advance counter
                V1loc_meridianRow = V1loc_meridianRow+1; % advance counter
            elseif any(~cellfun(@isempty,V1loc_meridian(V1loc_meridianRow,:))) % check that data is not being overwritten
                disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,3} ' at row ' num2str(V1loc_meridianRow) ' of V1_fixloc']) % message
                V1loc_meridianRow = V1loc_meridianRow+1; % advance counter
            end

            % store stats about eye tracking data
            V1loc_meridian{V1loc_meridianRow,1} = subjects{iS}; % store subject name
            V1loc_meridian{V1loc_meridianRow,2} = num2str(setN); % store set#
            V1loc_meridian{V1loc_meridianRow,3} = FDrms(iSet);
            V1loc_meridian{V1loc_meridianRow,4} = FDsd(iSet);
            V1loc_meridian{V1loc_meridianRow,5} = FDmax(iSet);

            % remove corresponding row from SubjPRTlist
            SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

        elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,3},V1locStr)==1 % is this a V1 localizer?

            % checks to advance row counter and avoid overwriting
            if all(~cellfun(@isempty,V1loc(V1locRow,:))) % check to advance counter
                V1locRow = V1locRow+1; % advance counter
            elseif any(~cellfun(@isempty,V1loc(V1locRow,:))) % check that data is not being overwritten
                disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,3} ' at row ' num2str(V1locRow) ' of V1loc']) % message
                V1locRow = V1locRow+1; % advance counter
            end

            % store stats about eye tracking data
            V1loc{V1locRow,1} = subjects{iS}; % store subject name
            V1loc{V1locRow,2} = num2str(setN); % store set#
            V1loc{V1locRow,3} = FDrms(iSet);
            V1loc{V1locRow,4} = FDsd(iSet);
            V1loc{V1locRow,5} = FDmax(iSet);

            % remove corresponding row from SubjPRTlist
            SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

        elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,3},supStr)==1 % is this a suppression scan?

            % checks to advance row counter and avoid overwriting
            if all(~cellfun(@isempty,sup(supRow,:))) % check to advance counter
                supRow = supRow+1; % advance counter
            elseif any(~cellfun(@isempty,sup(supRow,:))) % check that data is not being overwritten
                disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,3} ' at row ' num2str(supRow) ' of sup']) % message
                supRow = supRow+1; % advance counter
            end

            % store stats about eye tracking data
            sup{supRow,1} = subjects{iS}; % store subject name
            sup{supRow,2} = num2str(setN); % store set#
            sup{supRow,3} = FDrms(iSet);
            sup{supRow,4} = FDsd(iSet);
            sup{supRow,5} = FDmax(iSet);

            % remove corresponding row from SubjPRTlist
            SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

        elseif AK_whichPattern(SubjPRTlist{SubjPRTlistRow,3},sumStr)==1 % is this a summation scan?

            % checks to advance row counter and avoid overwriting
            if all(~cellfun(@isempty,sum(sumRow,:))) % check to advance counter
                sumRow = sumRow+1; % advance counter
            elseif any(~cellfun(@isempty,sum(sumRow,:))) % check that data is not being overwritten
                disp(['Avoided overwriting for set with associated prt file: ' SubjPRTlist{SubjPRTlistRow,3} ' at row ' num2str(sumRow) ' of sum']) % message
                sumRow = sumRow+1; % advance counter
            end

            % store stats about eye tracking data
            sum{sumRow,1} = subjects{iS}; % store subject name
            sum{sumRow,2} = num2str(setN); % store set#
            sum{sumRow,3} = FDrms(iSet);
            sum{sumRow,4} = FDsd(iSet);
            sum{sumRow,5} = FDmax(iSet);

            % remove corresponding row from SubjPRTlist
            SubjPRTlist(SubjPRTlistRow,:) = [];% remove row from SubjPRTlist

        end    
    end
    if ~isempty(SubjPRTlist)
        disp('***Some unmatched set#s remaining in prtList:'); % message
        SubjPRTlist
    end
end

% trim empty rows from end of each cell array
MTloc(all(cellfun(@isempty,MTloc(:,1:2)),2),:) = [];
V1loc(all(cellfun(@isempty,V1loc(:,1:2)),2),:) = [];
V1loc_meridian(all(cellfun(@isempty,V1loc_meridian(:,1:2)),2),:) = [];
sup(all(cellfun(@isempty,sup(:,1:2)),2),:) = [];
sum(all(cellfun(@isempty,sum(:,1:2)),2),:) = [];


%% write cell arrays to xls file spreadsheets

success = zeros(5,1); % preallocate
xlsFilename = fullfile(save_dir,'fMRI_HeadMotion_Data_Tables_SupSum.xlsx'); % name file for only subjects with good data

success(1) = xlswrite(xlsFilename,MTloc,'MT Loc');
success(2) = xlswrite(xlsFilename,V1loc,'V1 Loc');
success(3) = xlswrite(xlsFilename,V1loc_meridian,'V1_meridian Loc');
success(4) = xlswrite(xlsFilename,sup,'Suppression');
success(5) = xlswrite(xlsFilename,sum,'Summation');
