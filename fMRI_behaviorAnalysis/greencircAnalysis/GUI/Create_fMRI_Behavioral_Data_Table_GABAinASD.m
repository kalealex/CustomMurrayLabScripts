% Create_fMRI_Behavioral_Data_Table_GABAinASD
% Murray_Lab_2016
% Created by AMK on 7/7/16

%% designate directory, subject, and logfile indentifying information

top_dir = 'L:\MurrayLab\ASD\Data';
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab';
save_dir = 'L:\MurrayLab\DataTablesForGUI';
subjects = dir(fullfile(top_dir,'G*')); % look for folders matching subject code format
subjects = {subjects.name}; % cat names
subjects = [subjects {'KG122'}]; % add kiddo
% list current list of subjects here; last updated 7/7/16
% subjects = {'G101','G102','G103','G104','G105','G106','G107','G109','G110','G111','G112','G307','G310','G311','G312','G313','G314','G315','G316','G317','G318','G319','G320','G322','G327'};

% strings for logfile recognition
experimentStr = {'MTlocalizer','V1localizer','V1localizer_fix','contrast','suppression','summation'};

%% designate table structure and field names

% name columns of behavioral data tables
general_columns = {'subject','session','fMRI set#','logfile scan#'};
moving_columns = {'moving: hit rate','moving: false alarm rate'};
static_columns = {'static: hit rate','static: false alarm rate'};
big_columns = {'big: hit rate','big: false alarm rate'};
small_columns = {'small: hit rate','small: false alarm rate'};
fix_columns = {'fix: hit rate','fix: false alarm rate'};
low_columns = {'low: hit rate','low: false alarm rate'};
high_columns = {'high: hit rate','high: false alarm rate'};

% create cell arrays to store data and save into spreadsheets
MTloc = cell(3*length(subjects)+50,length(general_columns)+length(moving_columns)+length(static_columns)); % 50 rows for buffer
V1loc = cell(3*length(subjects)+50,length(general_columns)+length(big_columns)+length(small_columns)); % 50 rows for buffer
V1_fixloc = cell(3*length(subjects)+50,length(general_columns)+length(fix_columns)+length(small_columns)); % 50 rows for buffer
contrast = cell(6*length(subjects)+100,length(general_columns)+length(fix_columns)+length(low_columns)+length(high_columns)); % 100 rows for buffer
sup = cell(6*length(subjects)+100,length(general_columns)+length(small_columns)+length(big_columns)); % 100 rows for buffer; sup and sum have the same condition names but separate cell arrays
sum = cell(6*length(subjects)+100,length(general_columns)+length(small_columns)+length(big_columns)); % 100 rows for buffer

% add column names to cell arrays
MTloc(1,:) = [general_columns,moving_columns,static_columns];
V1loc(1,:) = [general_columns,big_columns,small_columns];
V1_fixloc(1,:) = [general_columns,fix_columns,small_columns];
contrast(1,:) = [general_columns,fix_columns,low_columns,high_columns];
sup(1,:) = [general_columns,small_columns,big_columns];
sum(1,:) = [general_columns,small_columns,big_columns];

% create counters for rows of each cell array
MTlocRow = 2;
V1locRow = 2;
V1_fixlocRow = 2;
contrastRow = 2;
supRow = 2;
sumRow = 2;

%% run function to parse behavioral data

% behavioral_data_filename = 'fMRI_Behavioral_Data.mat';
% behavioral_data_dir = fullfile(home_dir,behavioral_data_filename);
% if ~exist(behavioral_data_dir,'file') % look for saved data
    data = AK_GABAinASD_BehaviorGreenCirc(subjects,top_dir,home_dir);
%     save(behavioral_data_dir,'data');
% else
%     disp(['loading ' behavioral_data_dir]) % message
%     load(behavioral_data_dir,'data');
% end

%% get list of assicated prt filenames for each set

prtList = AK_GABAinASD_getAllAssociatedPRTnamesBySetSessionSubject(subjects,0);

%% pull data from data structure into appropriat spreadsheets

for iS = 1:length(data) % cycle through subjects 
    for iF = 1:length(data(iS).fMRIsession) % cycle through fMRI sessions
        if ~isempty(data(iS).fMRIsession{iF})
            % create prtList for the current subject and session
            clear SubjPRTlist
            SubjPRTlist = prtList((AK_findStrMatch(prtList(:,1),subjects{iS}) & AK_findStrMatch(prtList(:,2),data(iS).fMRIsession{iF})),:);

            % only bother trying to match up sets to .log files if the .prt file links hav been generated
            if ~isempty(SubjPRTlist) && ~all(cellfun(@isempty,SubjPRTlist(:,4)))
                for iL = 1:length(data(iS).logFile(iF,:)) % cycle through logfiles
                    clear patternN scanNindex BADindex scanNfindStr scanN
                    scanNfindStr = [subjects{iS} '_']; % string to be used to find scan# in logfile name

                    if any(AK_whichPattern(data(iS).logFile{iF,iL},experimentStr)) % is it a localizer scan?

                        %determine type of localizer scan
                        patternN = AK_whichPattern(data(iS).logFile{iF,iL},experimentStr);
                        if patternN(3) == 1 % is V1_fix? (must check this before V1 due to overlap in string name)

                            % find scanNfindStr; otherwise skip this logfile
                            % also skip if logfile doesn't have a complete pulse count
                            % also skip if 'BAD'
                            scanNindex = strfind(data(iS).logFile{iF,iL},scanNfindStr);
                            BADindex = strfind(data(iS).logFile{iF,iL},'BAD');
                            if ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)>=80 && isempty(BADindex) % should have 80 pulses

                                %find and store scan# as variable scanN
                                if ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1))) % check for two digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1));
                                elseif ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)))) % check for one digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)));
                                end

                                % checks to advance row counter and avoid overwriting
                                if all(~cellfun(@isempty,V1_fixloc(V1_fixlocRow,5:8))) % check to advance counter
%                                     loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
%                                     loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                                    V1_fixlocRow = V1_fixlocRow+1; % advance counter
                                elseif any(~cellfun(@isempty,V1_fixloc(V1_fixlocRow,5:8))) % check that data is not being overwritten
                                    disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(V1_fixlocRow) ' of V1_fixloc']) % message
%                                     loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
%                                     loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                                    V1_fixlocRow = V1_fixlocRow+1; % advance counter
                                end

                                % store hit and false alarm rate data
                                V1_fixloc{V1_fixlocRow,1} = data(iS).subject; % store subject name
                                V1_fixloc{V1_fixlocRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                                V1_fixloc{V1_fixlocRow,4} = num2str(scanN); % store scan#
                                V1_fixloc(V1_fixlocRow,5) = table2cell(data(iS).accuracy{iF,iL}('surr','hitRate')); % store hit rate for surr condition
                                V1_fixloc(V1_fixlocRow,6) = table2cell(data(iS).accuracy{iF,iL}('surr','falseAlarmRate')); % store false alarm rate for surr condition
                                V1_fixloc(V1_fixlocRow,7) = table2cell(data(iS).accuracy{iF,iL}('tar','hitRate')); % store hit rate for tar condition
                                V1_fixloc(V1_fixlocRow,8) = table2cell(data(iS).accuracy{iF,iL}('tar','falseAlarmRate')); % store false alarm rate for tar condition

                                % store appropriate set number and remove corresponding row from SubjPRTlist
                                clear SubjPRTlistScanIndex
                                SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{3},1),1); % create index to find set#
                                V1_fixloc{V1_fixlocRow,3} = ['set' num2str(SubjPRTlist{SubjPRTlistScanIndex,3})]; % store set#
                                SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist

                            elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<80 % can find scanNindex, but too few pulses
                                % find and comare the number of remaining sets and logfiles of the appropriate type
                                clear nSets nLogfiles
                                nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{3},1)));
                                nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),experimentStr{3})));

                                % determine whether or not the current set of this type needs to be removed from the list
                                if nSets >= nLogfiles 
                                    % if so, remove appropriate row from SubjPRTlist
                                    clear SubjPRTlistScanIndex
                                    SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{3},1),1); % create index to find set#
                                    SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist
                                end    
                            end
                        elseif patternN(2) == 1 % is V1?

                            % find scanNfindStr; otherwise skip this logfile
                            % also skip if logfile doesn't have a complete pulse count
                            % also skip if 'BAD'
                            scanNindex = strfind(data(iS).logFile{iF,iL},scanNfindStr);
                            BADindex = strfind(data(iS).logFile{iF,iL},'BAD');
                            if ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)>=80 && isempty(BADindex) % should have 80 pulses

                                %find and store scan# as variable scanN
                                if ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1))) % check for two digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1));
                                elseif ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)))) % check for one digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)));
                                end

                                % checks to advance row counter and avoid overwriting
                                if all(~cellfun(@isempty,V1loc(V1locRow,5:8))) % check to advance counter
%                                     loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
%                                     loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                                    V1locRow = V1locRow+1; % advance counter
                                elseif any(~cellfun(@isempty,V1loc(V1locRow,5:8))) % check that data is not being overwritten
                                    disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(V1locRow) ' of V1loc']) % message
%                                     loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
%                                     loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                                    V1locRow = V1locRow+1; % advance counter
                                end

                                % store hit and false alarm rate data
                                V1loc{V1locRow,1} = data(iS).subject; % store subject name
                                V1loc{V1locRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                                V1loc{V1locRow,4} = num2str(scanN); % store scan#
                                V1loc(V1locRow,5) = table2cell(data(iS).accuracy{iF,iL}('surr','hitRate')); % store hit rate for surr condition
                                V1loc(V1locRow,6) = table2cell(data(iS).accuracy{iF,iL}('surr','falseAlarmRate')); % store false alarm rate for surr condition
                                V1loc(V1locRow,7) = table2cell(data(iS).accuracy{iF,iL}('tar','hitRate')); % store hit rate for tar condition
                                V1loc(V1locRow,8) = table2cell(data(iS).accuracy{iF,iL}('tar','falseAlarmRate')); % store false alarm rate for tar condition

                                % store appropriate set number and remove corresponding row from SubjPRTlist
                                clear SubjPRTlistScanIndex
                                SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{2},1),1); % create index to find set#
                                V1loc{V1locRow,3} = ['set' num2str(SubjPRTlist{SubjPRTlistScanIndex,3})]; % store set#
                                SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist

                            elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<80 % can find scanNindex, but too few pulses
                                % find and comare the number of remaining sets and logfiles of the appropriate type
                                clear nSets nLogfiles
                                nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{2},1)));
                                nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),experimentStr{2})));

                                % determine whether or not the current set of this type needs to be removed from the list
                                if nSets >= nLogfiles 
                                    % if so, remove appropriate row from SubjPRTlist
                                    clear SubjPRTlistScanIndex
                                    SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{2},1),1); % create index to find set#
                                    SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist
                                end        
                            end
                        elseif patternN(1) == 1 % is MT?

                            % find scanNfindStr; otherwise skip this logfile
                            % also skip if logfile doesn't have a complete pulse count
                            % also skip if 'BAD'
                            scanNindex = strfind(data(iS).logFile{iF,iL},scanNfindStr);
                            BADindex = strfind(data(iS).logFile{iF,iL},'BAD');
                            if ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)>=125 && isempty(BADindex) % should have 125 pulses

                                %find and store scan# as variable scanN
                                if ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1))) % check for two digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1));
                                elseif ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)))) % check for one digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)));
                                end

                                % checks to advance row counter and avoid overwriting
                                if all(~cellfun(@isempty,MTloc(MTlocRow,5:8))) % check to advance counter
%                                     loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
%                                     loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                                    MTlocRow = MTlocRow+1; % advance counter
                                elseif any(~cellfun(@isempty,MTloc(MTlocRow,5:8))) % check that data is not being overwritten
                                    disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(MTlocRow) ' of MTloc']) % message
%                                     loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
%                                     loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                                    MTlocRow = MTlocRow+1; % advance counter
                                end

                                % store hit and false alarm rate data
                                MTloc{MTlocRow,1} = data(iS).subject; % store subject name
                                MTloc{MTlocRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                                MTloc{MTlocRow,4} = num2str(scanN); % store scan#
                                MTloc(MTlocRow,5) = table2cell(data(iS).accuracy{iF,iL}('moving','hitRate')); % store hit rate for moving condition
                                MTloc(MTlocRow,6) = table2cell(data(iS).accuracy{iF,iL}('moving','falseAlarmRate')); % store false alarm rate for moving condition
                                MTloc(MTlocRow,7) = table2cell(data(iS).accuracy{iF,iL}('static','hitRate')); % store hit rate for static condition
                                MTloc(MTlocRow,8) = table2cell(data(iS).accuracy{iF,iL}('static','falseAlarmRate')); % store false alarm rate for static condition

                                % store appropriate set number and remove corresponding row from SubjPRTlist
                                clear SubjPRTlistScanIndex
                                SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{1},1),1); % create index to find set#
                                MTloc{MTlocRow,3} = ['set' num2str(SubjPRTlist{SubjPRTlistScanIndex,3})]; % store set#
                                SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist

                            elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<125 % can find scanNindex, but too few pulses
                                % find and comare the number of remaining sets and logfiles of the appropriate type
                                clear nSets nLogfiles
                                nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{1},1)));
                                nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),experimentStr{1})));

                                % determine whether or not the current set of this type needs to be removed from the list
                                if nSets >= nLogfiles 
                                    % if so, remove appropriate row from SubjPRTlist
                                    clear SubjPRTlistScanIndex
                                    SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{1},1),1); % create index to find set#
                                    SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist
                                end
                            end
                        elseif patternN(4) == 1 % is it a contrast scan?

                            % find scanNfindStr; otherwise skip this logfile
                            % also skip if logfile doesn't have a complete pulse count
                            % also skip if 'BAD'
                            scanNindex = strfind(data(iS).logFile{iF,iL},scanNfindStr);
                            BADindex = strfind(data(iS).logFile{iF,iL},'BAD');
                            if ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)>=125 && isempty(BADindex) % should have 125 pulses

                                %find and store scan# as variable scanN
                                if ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1))) % check for two digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1));
                                elseif ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)))) % check for one digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)));
                                end

                                % store hit and false alarm rate data
                                contrast{contrastRow,1} = data(iS).subject; % store subject name
                                contrast{contrastRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                                contrast{contrastRow,4} = num2str(scanN); % store scan#
                                contrast(contrastRow,5) = table2cell(data(iS).accuracy{iF,iL}('fix','hitRate')); % store hit rate for fix condition
                                contrast(contrastRow,6) = table2cell(data(iS).accuracy{iF,iL}('fix','falseAlarmRate')); % store false alarm rate for fix condition
                                contrast(contrastRow,7) = table2cell(data(iS).accuracy{iF,iL}('low','hitRate')); % store hit rate for low condition
                                contrast(contrastRow,8) = table2cell(data(iS).accuracy{iF,iL}('low','falseAlarmRate')); % store false alarm rate for low condition
                                contrast(contrastRow,9) = table2cell(data(iS).accuracy{iF,iL}('hi','hitRate')); % store hit rate for hi condition
                                contrast(contrastRow,10) = table2cell(data(iS).accuracy{iF,iL}('hi','falseAlarmRate')); % store false alarm rate for hi condition

                                % store appropriate set number and remove corresponding row from SubjPRTlist
                                clear SubjPRTlistScanIndex
                                SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{4},1),1); % create index to find set#
                                contrast{contrastRow,3} = ['set' num2str(SubjPRTlist{SubjPRTlistScanIndex,3})]; % store next set# of scan type
                                SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist

                                contrastRow = contrastRow+1; % advance counter

                            elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<125 % can find scanNindex, but too few pulses
                                % find and comare the number of remaining sets and logfiles of the appropriate type
                                clear nSets nLogfiles
                                nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{4},1)));
                                nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),experimentStr{4})));

                                % determine whether or not the current set of this type needs to be removed from the list
                                if nSets>=nLogfiles 
                                    % if so, remove appropriate row from SubjPRTlist
                                    clear SubjPRTlistScanIndex
                                    SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{4},1)); % create index to find set#
                                    SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                                end       
                            end
                        elseif patternN(5) == 1 % is it a suppression scan?

                            % find scanNfindStr; otherwise skip this logfile
                            % also skip if logfile doesn't have a complete pulse count
                            % also skip if 'BAD'
                            scanNindex = strfind(data(iS).logFile{iF,iL},scanNfindStr);
                            BADindex = strfind(data(iS).logFile{iF,iL},'BAD');
                            if ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)>=125 && isempty(BADindex) % should have 125 pulses

                                %find and store scan# as variable scanN
                                if ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1))) % check for two digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1));
                                elseif ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)))) % check for one digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)));
                                end

                                % checks to advance row counter and avoid overwriting
                                if all(~cellfun(@isempty,sup(supRow,5:8))) % check to advance counter
%                                     supsum{supsumRow,3} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,3},'UniformOutput',false),','); % convert set# to string
%                                     supsum{supsumRow,4} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,4},'UniformOutput',false),','); % convert scan# to string
                                    supRow = supRow+1; % advance counter
                                elseif any(~cellfun(@isempty,sup(supRow,5:8))) % check that data is not being overwritten
                                    disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(supRow) ' of sup']) % message
%                                     supsum{supsumRow,3} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,3},'UniformOutput',false),','); % convert set# to string
%                                     supsum{supsumRow,4} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,4},'UniformOutput',false),','); % convert scan# to string
                                    supRow = supRow+1; % advance counter
                                end

                                % store hit and false alarm rate data
                                sup{supRow,1} = data(iS).subject; % store subject name
                                sup{supRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                                sup{supRow,4} = num2str(scanN); % store scan#
                                sup(supRow,5) = table2cell(data(iS).accuracy{iF,iL}('small','hitRate')); % store hit rate for small condition
                                sup(supRow,6) = table2cell(data(iS).accuracy{iF,iL}('small','falseAlarmRate')); % store false alarm rate for small condition
                                sup(supRow,7) = table2cell(data(iS).accuracy{iF,iL}('big','hitRate')); % store hit rate for big condition
                                sup(supRow,8) = table2cell(data(iS).accuracy{iF,iL}('big','falseAlarmRate')); % store false alarm rate for big condition

                                % store appropriate set number and remove corresponding row from SubjPRTlist
                                clear SubjPRTlistScanIndex
                                SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{5},1),1); % create index to find set#
                                sup{supRow,3} = ['set' num2str(SubjPRTlist{SubjPRTlistScanIndex,3})]; % store set#
                                SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist

                            elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<125 % can find scanNindex, but too few pulses
                                % find and comare the number of remaining sets and logfiles of the appropriate type
                                clear nSets nLogfiles
                                nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{5},1)));
                                nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),experimentStr{5})));

                                % determine whether or not the current set of this type needs to be removed from the list
                                if nSets >= nLogfiles 
                                    % if so, remove appropriate row from SubjPRTlist
                                    clear SubjPRTlistScanIndex
                                    SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{5},1),1); % create index to find set#
                                    SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist
                                end
                            end
                        elseif patternN(6) == 1 % is suppression?

                            % find scanNfindStr; otherwise skip this logfile
                            % also skip if logfile doesn't have a complete pulse count
                            % also skip if 'BAD'
                            scanNindex = strfind(data(iS).logFile{iF,iL},scanNfindStr);
                            BADindex = strfind(data(iS).logFile{iF,iL},'BAD');
                            if ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)>=125 && isempty(BADindex) % should have 125 pulses

                                %find and store scan# as variable scanN
                                if ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1))) % check for two digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr):scanNindex+length(scanNfindStr)+1));
                                elseif ~isnan(str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)))) % check for one digit scan#
                                    scanN = str2double(data(iS).logFile{iF,iL}(scanNindex+length(scanNfindStr)));
                                end

                                % checks to advance row counter and avoid overwriting
                                if all(~cellfun(@isempty,sum(sumRow,5:8))) % check to advance counter
%                                     supsum{supsumRow,3} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,3},'UniformOutput',false),','); % convert set# to string
%                                     supsum{supsumRow,4} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,4},'UniformOutput',false),','); % convert scan# to string
                                    sumRow = sumRow+1; % advance counter
                                elseif any(~cellfun(@isempty,sum(sumRow,5:8))) % check that data is not being overwritten
                                    disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(sumRow) ' of sum']) % message
%                                     supsum{supsumRow,3} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,3},'UniformOutput',false),','); % convert set# to string
%                                     supsum{supsumRow,4} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,4},'UniformOutput',false),','); % convert scan# to string
                                    sumRow = sumRow+1; % advance counter
                                end

                                % store hit and false alarm rate data
                                sum{sumRow,1} = data(iS).subject; % store subject name
                                sum{sumRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                                sum{sumRow,4} = num2str(scanN); % store scan#
                                sum(sumRow,5) = table2cell(data(iS).accuracy{iF,iL}('small','hitRate')); % store hit rate for small condition
                                sum(sumRow,6) = table2cell(data(iS).accuracy{iF,iL}('small','falseAlarmRate')); % store false alarm rate for small condition
                                sum(sumRow,7) = table2cell(data(iS).accuracy{iF,iL}('big','hitRate')); % store hit rate for big condition
                                sum(sumRow,8) = table2cell(data(iS).accuracy{iF,iL}('big','falseAlarmRate')); % store false alarm rate for big condition

                                % store appropriate set number and remove corresponding row from SubjPRTlist
                                clear SubjPRTlistScanIndex
                                SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{6},1),1); % create index to find set#
                                sum{sumRow,3} = ['set' num2str(SubjPRTlist{SubjPRTlistScanIndex,3})]; % store set#
                                SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist

                            elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<125 % can find scanNindex, but too few pulses
                                % find and comare the number of remaining sets and logfiles of the appropriate type
                                clear nSets nLogfiles
                                nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{6},1)));
                                nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),experimentStr{6})));

                                % determine whether or not the current set of this type needs to be removed from the list
                                if nSets >= nLogfiles 
                                    % if so, remove appropriate row from SubjPRTlist
                                    clear SubjPRTlistScanIndex
                                    SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),experimentStr{6},1),1); % create index to find set#
                                    SubjPRTlist(SubjPRTlistScanIndex,:) = [];% remove row from SubjPRTlist
                                end    
                            end
                        end
                    end
                end
            end
        end
        if ~isempty(SubjPRTlist)
            disp('***Some unmatched set#s remaining in prtList:'); % message
            SubjPRTlist
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

%% write cell arrays to xls file spreadsheets

success = zeros(6,1); % preallocate
xlsFilename = fullfile(save_dir,'fMRI_Behavioral_Data_Tables.xlsx'); % name file

success(1) = xlswrite(xlsFilename,MTloc,'MT_localizer');
success(2) = xlswrite(xlsFilename,V1loc,'V1_localizer');
success(3) = xlswrite(xlsFilename,V1_fixloc,'V1_fix_localizer');
success(4) = xlswrite(xlsFilename,contrast,'Contrast');
success(5) = xlswrite(xlsFilename,sup,'Suppression');
success(6) = xlswrite(xlsFilename,sum,'Summation');


