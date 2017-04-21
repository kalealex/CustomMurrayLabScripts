% Create_fMRI_Behavioral_Data_Table_GABAinASD
% Murray_Lab_2016
% Created by AMK on 7/7/16

%% designate directory, subject, and logfile indentifying information

top_dir = 'L:\MurrayLab\ASD\Data';
home_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab';
% list current list of subjects here; last updated 7/7/16
subjects = {'G101','G102','G103','G104','G105','G106','G107','G109','G110','G111','G112','G307','G310','G311','G312','G313','G314','G315','G316','G317','G318','G319','G320','G322','G327'};

% strings for logfile recognition
locStr = {'MTlocalizer','V1localizer','V1localizer_fix'};
contrastStr = {'contrast'};
supsumStr = {'suppression','summation'};

%% designate table structure and field names

% name columns of behavioral data tables
general_columns = {'subject','session','fMRI set#','logfile scan#'};
loc_columns = {'MT moving: hit rate','MT moving: false alarm rate','MT static: hit rate','MT static: false alarm rate','V1 big: hit rate','V1 big: false alarm rate','V1 small: hit rate','V1 small: false alarm rate','V1_fix fix: hit rate','V1_fix fix: false alarm rate','V1_fix small: hit rate','V1_fix small: false alarm rate'};
contrast_columns = {'fix: hit rate','fix: false alarm rate','low: hit rate','low: false alarm rate','high: hit rate','high: false alarm rate'};
supsum_columns = {'sup small: hit rate','sup small: false alarm rate','sup big: hit rate','sup big: false alarm rate','sum small: hit rate','sum small: false alarm rate','sum big: hit rate','sum big: false alarm rate'};

% create cell arrays to store data and save into spreadsheets
loc = cell(2*length(subjects)+50,length(general_columns)+length(loc_columns)); % 50 rows for buffer
contrast = cell(4*length(subjects)+100,length(general_columns)+length(contrast_columns)); % 100 rows for buffer
supsum = cell(4*length(subjects)+100,length(general_columns)+length(supsum_columns)); % 100 rows for buffer

% add column names to cell arrays
loc(1,:) = [general_columns,loc_columns];
contrast(1,:) = [general_columns,contrast_columns];
supsum(1,:) = [general_columns,supsum_columns];

% set up set# & scan# index order that matches table organization
loc(2:end,3:4) = {zeros(1,3)}; % one place for each type of scan (MT,V1,V1_fix)
supsum(2:end,3:4) = {zeros(1,2)}; % one place for each type of scan (sup,sum)


% create counters for rows of each cell array
locRow = 2;
contrastRow = 2;
supsumRow = 2;

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
        % create prtList for the current subject and session
        clear SubjPRTlist
        SubjPRTlist = prtList((AK_findStrMatch(prtList(:,1),subjects{iS}) & AK_findStrMatch(prtList(:,2),data(iS).fMRIsession{iF})),:); 
        for iL = 1:length(data(iS).logFile(iF,:)) % cycle through logfiles
            clear patternN scanNindex BADindex scanNfindStr scanN
            scanNfindStr = [subjects{iS} '_']; % string to be used to find scan# in logfile name
            
            if any(AK_whichPattern(data(iS).logFile{iF,iL},locStr)) % is it a localizer scan?
                
                %determine type of localizer scan
                patternN = AK_whichPattern(data(iS).logFile{iF,iL},locStr);
                if patternN(3)==1 % is V1_fix? (must check this before V1 due to overlap in string name)
                    
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
                        if all(~cellfun(@isempty,loc(locRow,5:16))) % check to advance counter
                            loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
                            loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                            locRow = locRow+1; % advance counter
                        elseif any(~cellfun(@isempty,loc(locRow,13:16))) % check that data is not being overwritten
                            disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(locRow) ' of loc']) % message
                            loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
                            loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                            locRow = locRow+1; % advance counter
                        end
                        
                        % store hit and false alarm rate data
                        loc{locRow,1} = data(iS).subject; % store subject name
                        loc{locRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                        loc{locRow,4}(3) = scanN; % store scan# at 3rd position (V1_fix)
                        loc(locRow,13) = table2cell(data(iS).accuracy{iF,iL}('surr','hitRate')); % store hit rate for surr condition
                        loc(locRow,14) = table2cell(data(iS).accuracy{iF,iL}('surr','falseAlarmRate')); % store false alarm rate for surr condition
                        loc(locRow,15) = table2cell(data(iS).accuracy{iF,iL}('tar','hitRate')); % store hit rate for tar condition
                        loc(locRow,16) = table2cell(data(iS).accuracy{iF,iL}('tar','falseAlarmRate')); % store false alarm rate for tar condition
                        
                        % store appropriate set number and remove corresponding row from SubjPRTlist
                        clear SubjPRTlistScanIndex
                        SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),locStr{3},1)); % create index to find set#
                        loc{locRow,3}(3) = SubjPRTlist{SubjPRTlistScanIndex(1),3}; % store next set# of scan type at 3rd position (V1_fix)
                        SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                        
                    elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<80 % can find scanNindex, but too few pulses
                        % find and comare the number of remaining sets and logfiles of the appropriate type
                        clear nSets nLogfiles
                        nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),locStr{3},1)));
                        nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),locStr{3})));
                        
                        % determine whether or not the current set of this type needs to be removed from the list
                        if nSets>=nLogfiles 
                            % if so, remove appropriate row from SubjPRTlist
                            clear SubjPRTlistScanIndex
                            SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),locStr{3},1)); % create index to find set#
                            SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                        end    
                    end
                elseif patternN(2)==1 % is V1?
                    
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
                        if all(~cellfun(@isempty,loc(locRow,5:16))) % check to advance counter
                            loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
                            loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                            locRow = locRow+1; % advance counter
                        elseif any(~cellfun(@isempty,loc(locRow,9:12))) % check that data is not being overwritten
                            disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(locRow) ' of loc']) % message
                            loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
                            loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                            locRow = locRow+1; % advance counter
                        end
                        
                        % store hit and false alarm rate data
                        loc{locRow,1} = data(iS).subject; % store subject name
                        loc{locRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                        loc{locRow,4}(2) = scanN; % store scan# at 2nd position (V1)
                        loc(locRow,9) = table2cell(data(iS).accuracy{iF,iL}('surr','hitRate')); % store hit rate for surr condition
                        loc(locRow,10) = table2cell(data(iS).accuracy{iF,iL}('surr','falseAlarmRate')); % store false alarm rate for surr condition
                        loc(locRow,11) = table2cell(data(iS).accuracy{iF,iL}('tar','hitRate')); % store hit rate for tar condition
                        loc(locRow,12) = table2cell(data(iS).accuracy{iF,iL}('tar','falseAlarmRate')); % store false alarm rate for tar condition
                        
                        % store appropriate set number and remove corresponding row from SubjPRTlist
                        clear SubjPRTlistScanIndex
                        SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),locStr{2},1)); % create index to find set#
                        loc{locRow,3}(2) = SubjPRTlist{SubjPRTlistScanIndex(1),3}; % store next set# of scan type at 2nd position (V1)
                        SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist

                    elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<80 % can find scanNindex, but too few pulses
                        % find and comare the number of remaining sets and logfiles of the appropriate type
                        clear nSets nLogfiles
                        nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),locStr{2},1)));
                        nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),locStr{2})));
                        
                        % determine whether or not the current set of this type needs to be removed from the list
                        if nSets>=nLogfiles 
                            % if so, remove appropriate row from SubjPRTlist
                            clear SubjPRTlistScanIndex
                            SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),locStr{2},1)); % create index to find set#
                            SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                        end        
                    end
                elseif patternN(1)==1 % is MT?
                    
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
                        if all(~cellfun(@isempty,loc(locRow,5:16))) % check to advance counter
                            loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
                            loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                            locRow = locRow+1; % advance counter
                        elseif any(~cellfun(@isempty,loc(locRow,5:8))) % check that data is not being overwritten
                            disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(locRow) ' of loc']) % message
                            loc{locRow,3} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,3},'UniformOutput',false),','); % convert set# to string
                            loc{locRow,4} = strjoin(arrayfun(@(x) num2str(x),loc{locRow,4},'UniformOutput',false),','); % convert scan# to string
                            locRow = locRow+1; % advance counter
                        end
                        
                        % store hit and false alarm rate data
                        loc{locRow,1} = data(iS).subject; % store subject name
                        loc{locRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                        loc{locRow,4}(1) = scanN; % store scan# at 1st position (MT)
                        loc(locRow,5) = table2cell(data(iS).accuracy{iF,iL}('moving','hitRate')); % store hit rate for moving condition
                        loc(locRow,6) = table2cell(data(iS).accuracy{iF,iL}('moving','falseAlarmRate')); % store false alarm rate for moving condition
                        loc(locRow,7) = table2cell(data(iS).accuracy{iF,iL}('static','hitRate')); % store hit rate for static condition
                        loc(locRow,8) = table2cell(data(iS).accuracy{iF,iL}('static','falseAlarmRate')); % store false alarm rate for static condition
                        
                        % store appropriate set number and remove corresponding row from SubjPRTlist
                        clear SubjPRTlistScanIndex
                        SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),locStr{1},1)); % create index to find set#
                        loc{locRow,3}(1) = SubjPRTlist{SubjPRTlistScanIndex(1),3}; % store next set# of scan type at 1st position (MT)
                        SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist

                    elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<125 % can find scanNindex, but too few pulses
                        % find and comare the number of remaining sets and logfiles of the appropriate type
                        clear nSets nLogfiles
                        nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),locStr{1},1)));
                        nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),locStr{1})));
                        
                        % determine whether or not the current set of this type needs to be removed from the list
                        if nSets>=nLogfiles 
                            % if so, remove appropriate row from SubjPRTlist
                            clear SubjPRTlistScanIndex
                            SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),locStr{1},1)); % create index to find set#
                            SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                        end            
                    end
                end
            elseif any(AK_whichPattern(data(iS).logFile{iF,iL},contrastStr)) % is it a contrast scan?
                
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
                    contrast{contrastRow,4} = scanN; % store scan#
                    contrast(contrastRow,5) = table2cell(data(iS).accuracy{iF,iL}('fix','hitRate')); % store hit rate for fix condition
                    contrast(contrastRow,6) = table2cell(data(iS).accuracy{iF,iL}('fix','falseAlarmRate')); % store false alarm rate for fix condition
                    contrast(contrastRow,7) = table2cell(data(iS).accuracy{iF,iL}('low','hitRate')); % store hit rate for low condition
                    contrast(contrastRow,8) = table2cell(data(iS).accuracy{iF,iL}('low','falseAlarmRate')); % store false alarm rate for low condition
                    contrast(contrastRow,9) = table2cell(data(iS).accuracy{iF,iL}('hi','hitRate')); % store hit rate for hi condition
                    contrast(contrastRow,10) = table2cell(data(iS).accuracy{iF,iL}('hi','falseAlarmRate')); % store false alarm rate for hi condition
                    
                    % store appropriate set number and remove corresponding row from SubjPRTlist
                    clear SubjPRTlistScanIndex
                    SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),contrastStr{1},1)); % create index to find set#
                    contrast{contrastRow,3} = SubjPRTlist{SubjPRTlistScanIndex(1),3}; % store next set# of scan type
                    SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                    
                    contrast{contrastRow,3} = num2str(contrast{contrastRow,3}); % convert set# to string
                    contrast{contrastRow,4} = num2str(contrast{contrastRow,4}); % convert scan# to string
                    contrastRow = contrastRow+1; % advance counter
                
                elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<125 % can find scanNindex, but too few pulses
                    % find and comare the number of remaining sets and logfiles of the appropriate type
                    clear nSets nLogfiles
                    nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),contrastStr{1},1)));
                    nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),contrastStr{1})));

                    % determine whether or not the current set of this type needs to be removed from the list
                    if nSets>=nLogfiles 
                        % if so, remove appropriate row from SubjPRTlist
                        clear SubjPRTlistScanIndex
                        SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),contrastStr{1},1)); % create index to find set#
                        SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                    end       
                end
            elseif any(AK_whichPattern(data(iS).logFile{iF,iL},supsumStr)) % is it a supsum scan?
                
                %determine type of supsum scan
                patternN = AK_whichPattern(data(iS).logFile{iF,iL},supsumStr);
                if patternN(2)==1 % is summation? 
                    
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
                        if all(~cellfun(@isempty,supsum(supsumRow,5:12))) % check to advance counter
                            supsum{supsumRow,3} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,3},'UniformOutput',false),','); % convert set# to string
                            supsum{supsumRow,4} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,4},'UniformOutput',false),','); % convert scan# to string
                            supsumRow = supsumRow+1; % advance counter
                        elseif any(~cellfun(@isempty,supsum(supsumRow,9:12))) % check that data is not being overwritten
                            disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(supsumRow) ' of supsum']) % message
                            supsum{supsumRow,3} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,3},'UniformOutput',false),','); % convert set# to string
                            supsum{supsumRow,4} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,4},'UniformOutput',false),','); % convert scan# to string
                            supsumRow = supsumRow+1; % advance counter
                        end
                        
                        % store hit and false alarm rate data
                        supsum{supsumRow,1} = data(iS).subject; % store subject name
                        supsum{supsumRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                        supsum{supsumRow,4}(2) = scanN; % store scan# at 2nd position (sum)
                        supsum(supsumRow,9) = table2cell(data(iS).accuracy{iF,iL}('small','hitRate')); % store hit rate for small condition
                        supsum(supsumRow,10) = table2cell(data(iS).accuracy{iF,iL}('small','falseAlarmRate')); % store false alarm rate for small condition
                        supsum(supsumRow,11) = table2cell(data(iS).accuracy{iF,iL}('big','hitRate')); % store hit rate for big condition
                        supsum(supsumRow,12) = table2cell(data(iS).accuracy{iF,iL}('big','falseAlarmRate')); % store false alarm rate for big condition
                        
                        % store appropriate set number and remove corresponding row from SubjPRTlist
                        clear SubjPRTlistScanIndex
                        SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),supsumStr{2},1)); % create index to find set#
                        supsum{supsumRow,3}(2) = SubjPRTlist{SubjPRTlistScanIndex(1),3}; % store next set# of scan type at 2nd position (sum)
                        SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                    
                    elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<125 % can find scanNindex, but too few pulses
                        % find and comare the number of remaining sets and logfiles of the appropriate type
                        clear nSets nLogfiles
                        nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),supsumStr{2},1)));
                        nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),supsumStr{2})));
                        
                        % determine whether or not the current set of this type needs to be removed from the list
                        if nSets>=nLogfiles 
                            % if so, remove appropriate row from SubjPRTlist
                            clear SubjPRTlistScanIndex
                            SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),supsumStr{2},1)); % create index to find set#
                            SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                        end
                    end
                elseif patternN(1)==1 % is suppression?
                    
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
                        if all(~cellfun(@isempty,supsum(supsumRow,5:12))) % check to advance counter
                            supsum{supsumRow,3} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,3},'UniformOutput',false),','); % convert set# to string
                            supsum{supsumRow,4} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,4},'UniformOutput',false),','); % convert scan# to string
                            supsumRow = supsumRow+1; % advance counter
                        elseif any(~cellfun(@isempty,supsum(supsumRow,5:8))) % check that data is not being overwritten
                            disp(['Avoided overwriting for logfile: ' data(iS).logFile{iF,iL} ' at row ' num2str(supsumRow) ' of supsum']) % message
                            supsum{supsumRow,3} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,3},'UniformOutput',false),','); % convert set# to string
                            supsum{supsumRow,4} = strjoin(arrayfun(@(x) num2str(x),supsum{supsumRow,4},'UniformOutput',false),','); % convert scan# to string
                            supsumRow = supsumRow+1; % advance counter
                        end
                        
                        % store hit and false alarm rate data
                        supsum{supsumRow,1} = data(iS).subject; % store subject name
                        supsum{supsumRow,2} = data(iS).fMRIsession{iF}; % store fMRI session name
                        supsum{supsumRow,4}(1) = scanN; % store scan# at 1st position (sup)
                        supsum(supsumRow,5) = table2cell(data(iS).accuracy{iF,iL}('small','hitRate')); % store hit rate for small condition
                        supsum(supsumRow,6) = table2cell(data(iS).accuracy{iF,iL}('small','falseAlarmRate')); % store false alarm rate for small condition
                        supsum(supsumRow,7) = table2cell(data(iS).accuracy{iF,iL}('big','hitRate')); % store hit rate for big condition
                        supsum(supsumRow,8) = table2cell(data(iS).accuracy{iF,iL}('big','falseAlarmRate')); % store false alarm rate for big condition
                        
                        % store appropriate set number and remove corresponding row from SubjPRTlist
                        clear SubjPRTlistScanIndex
                        SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),supsumStr{1},1)); % create index to find set#
                        supsum{supsumRow,3}(1) = SubjPRTlist{SubjPRTlistScanIndex(1),3}; % store next set# of scan type at 1st position (sup)
                        SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
                        
                    elseif ~isempty(scanNindex) && data(iS).pulseCount(iF,iL)<125 % can find scanNindex, but too few pulses
                        % find and comare the number of remaining sets and logfiles of the appropriate type
                        clear nSets nLogfiles
                        nSets = length(find(AK_findStrMatch(SubjPRTlist(:,4),supsumStr{1},1)));
                        nLogfiles = length(find(AK_findStrMatch(data(iS).logFile(iF,iL:find(~cellfun(@isempty,data(iS).logFile(iF,:)),1,'last')),supsumStr{1})));
                        
                        % determine whether or not the current set of this type needs to be removed from the list
                        if nSets>=nLogfiles 
                            % if so, remove appropriate row from SubjPRTlist
                            clear SubjPRTlistScanIndex
                            SubjPRTlistScanIndex = find(AK_findStrMatch(SubjPRTlist(:,4),supsumStr{1},1)); % create index to find set#
                            SubjPRTlist(SubjPRTlistScanIndex(1),:) = [];% remove row from SubjPRTlist
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
loc(all(cellfun(@isempty,loc(:,1:2)),2),:) = [];
contrast(all(cellfun(@isempty,contrast(:,1:2)),2),:) = [];
supsum(all(cellfun(@isempty,supsum(:,1:2)),2),:) = [];

%% write cell arrays to xls file spreadsheets

success = zeros(3,1); % preallocate
xlsFilename = fullfile(top_dir,'fMRI_Behavioral_Data_Tables.xlsx'); % name file

success(1) = xlswrite(xlsFilename,loc,'Localizers');
success(2) = xlswrite(xlsFilename,contrast,'Contrast');
success(3) = xlswrite(xlsFilename,supsum,'SupSum');


