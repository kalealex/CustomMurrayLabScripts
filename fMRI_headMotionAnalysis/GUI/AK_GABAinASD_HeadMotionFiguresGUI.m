function AK_GABAinASD_HeadMotionFiguresGUI( xlsfile )
%AK_GABAinASD_HeadMotionFiguresGUI is a guided user interface designed to 
%read in specifically formatted data from a .xls file and turn it into a 
%beeswarm plot. Specifically, the function creates a subplot for each
%statistic selected via the GUI so that each experiment selected receives 
%its own distribution spaced out along the x axis. Within each experiment 
%type, data are further sparated into groups (ASD and NT) so that each 
%group * experiment combination has a distribution plotted on a common axis.
%There are also options to filter data to appear in different colors within 
%each distribution based on subject ID or based on a criterion and some 
%other statistic (i.e., the filter framewise displacement SD > .8 would 
%highlight points where framewise displacement SD exceeds .8 in a different
%color than points which do not me that criterion)
%
% Run the script Create_fMRI_HeadMotion_Data_Table_GABAinASD.m in order to
% generate correctly formatted .xls files for the GUI to read. The existing
% file on the L drive should be named fMRI_HeadMotion_Data_Tables*.xlsx.
%   
%INPUT:
%   xlsfile [optional]: full file directory for .xls file where head motion data
%       tables are stored (string); defaults to UI selection from
%       directory: 'L:\MurrayLab\DataTablesForGUI'
%OUTPUT: (option to save output using 'Export' button)
%   exportData: a maxrix with 5 dimensions (subject, session, set, 
%       experiment, statistic), the content of which will match whatever 
%       data set is currently chosen for visualization; blank indices are 
%       filled with nans 
%   key: a structure with fields containing labels for each dimesion of
%       exportData (i.e., key.d1 contains a list of subjects whose data is in
%       exportData)


%%% Changes from EyeTracking Version of GUI:
% this version of the GUI does not read conditions
    % simplified reading of column headers in Spreadsheads: they no longer need to be split into conds and stats; just stats now
    % remome conditions from x-axes and data indexing in RePlot
        % everything that was indexed by condition is now indexed by sheet (experiment)
% this version also doesn't deal with spreadsheets for block-level statistics
    % changed hard coded sheet names in sheetsAll
    % removed blockSheet statements and indexing


%% check input 

if nargin < 1
    xlsfile_dir = 'L:\MurrayLab\DataTablesForGUI'; % set root dir (should be the same for any cpu where network drive is mapped to 'L:\')
    cd(xlsfile_dir);
    xlsfile_name = uigetfile('*.xlsx','Select a file from which to load data');
    xlsfile = fullfile(xlsfile_dir,xlsfile_name); 
end

addpath(genpath('L:\MurrayLab\ASD\Data\Analysis_scripts\AMK_Code'));

%% initialize display and set up GUI referencing callback functions

% set prior values for variables spanning multiple functions
sheetsCurrent = {''}; % to pass between callback functions
% condsCurrent = {''};
statsCurrent = {''};
subjectsCurrent = {''};
filterIdx = {[]};
filtstatSelected = '';
operatorSelected = '';
criterionSelected = [];
highlight_subjects = {''};
plotStruct = struct;

dataTables = {[]}; % to be defined in LoadData
% isBlockSheet = [];
% condsAll = {''};
% condsAllIdx = {[]};
% condsAllWhichSheetIdx = [];
statsAll = {''};
statsAllIdx = {[]};
subjectsAll = {''};

exportData = []; % optional output (needs to have scope of full function)

% create or find figure
fH = findobj(0,'tag','GABAinASD_EyeTracking_Fig');
if ~isempty(fH)
    figure(fH);
else
    InitializeFigureWindow();
end
    
%% nested functions
    
    function InitializeFigureWindow()
        % create figure
        fH = figure('tag','GABAinASD_EyeTracking_Fig', 'name','GABAinASD_EyeTracking_Fig','toolbar','figure','numbertitle','off');
        % get figure position and size
        figPosition = getpixelposition(fH);
        % set up UI tags
        uicontrol('tag','Experiments','string','Experiments','callback',@Spreadsheets_CB,'pos',[5 112 70 20]);
        uicontrol('tag','Statistics','string','Statistics','callback',@Statistics_CB,'pos',[5 89 60 20]);
%         uicontrol('tag','Conditions','string','Conditions','callback',@Conditions_CB,'pos',[5 66 60 20]);
        uicontrol('tag','Subjects','string','Subjects','callback',@Subjects_CB,'pos',[5 66 60 20]);
        uicontrol('tag','Filter','string','Filter Data','callback',@Filter_CB,'pos',[5 43 60 20]);
        uicontrol('tag','ExportButton','string','Export','callback',@ExportData_CB,'pos',[5 20 60 20]);
    end

    function ExportData_CB(~,~)
        % get user input for where to save file
        [fileName,pathName] = uiputfile;
        % save formatted data
        if fileName ~= 0 % did they enter a name and directory?
%             if length(size(exportData)) == 6;
%                 key.d1 = subjectsCurrent;
%                 key.d2 = {'fMRI1','fMRI2'};
%                 key.d3 = {'set1','set2','set3','set4','set5','set6','set7','set8','set9'};
%                 key.d4 = 'blocks';
%                 key.d5 = condsCurrent;
%                 key.d6 = statsCurrent;
%             else
%                 key.d1 = subjectsCurrent;
%                 key.d2 = {'fMRI1','fMRI2'};
%                 key.d3 = {'set1','set2','set3','set4','set5','set6','set7','set8','set9'};
%                 key.d4 = condsCurrent;
%                 key.d5 = statsCurrent;
%             end
            key.d1 = subjectsCurrent;
            key.d2 = {'fMRI1','fMRI2'};
            key.d3 = {'set1','set2','set3','set4','set5','set6','set7','set8','set9'};
            key.d4 = sheetsCurrent; % where sheets are experiments
            key.d5 = statsCurrent;
            save(fullfile(pathName,fileName),'exportData','key')
        end
    end
    
    function LoadData() 
        % load eye tracking data tables
        for iSheet = 1:length(sheetsCurrent)
            % read in data
            [~,~,dataTables{iSheet}] = xlsread(xlsfile,sheetsCurrent{iSheet}); % read sheet
            % is this a sheet for block statistics?
%             isBlockSheet(iSheet) = ~isempty(regexp(sheetsCurrent{iSheet},regexptranslate('wildcard','block*'),'ONCE'));
            % fill in empty cells with nans
            dataTables{iSheet}(cellfun(@isempty,dataTables{iSheet})) = {nan};
            % eliminate columns of all nans or mostly nans (these represent notes taken on the spreadsheet file)
            dataTables{iSheet}(:,all(cellfun(@(x) all(isnan(x)),dataTables{iSheet}),1) | sum(~cellfun(@(x) all(isnan(x)),dataTables{iSheet}),1) < 5 ) = [];
            
            % find subjects list
            tableSubjects{iSheet} = unique(dataTables{iSheet}(2:end,1));
            % find condition and stat names for current sheet
%             clear condStrIndex tableCondListFull tableStatListFull tableCondListUnique condIdx statIdx tableCondList
%             findCondStatStr = '*:*';
%             condStatStrIndex = find(~cellfun(@isempty,regexp(dataTables{iSheet}(1,:),regexptranslate('wildcard',findCondStatStr)))==1);
%             for iCS = 1:length(condStatStrIndex)
%                 clear tempIdx
%                 tempIdx = strfind(dataTables{iSheet}{1,condStatStrIndex(iCS)},findCondStatStr(2)); % idx for position of : in each string where : occurs
%                 tableCondListFull{iCS} = dataTables{iSheet}{1,condStatStrIndex(iCS)}(1:tempIdx-1); % list of all condition names in order of columns
%                 tableStatListFull{iCS} = dataTables{iSheet}{1,condStatStrIndex(iCS)}(tempIdx+2:end); % list of all stat names in order of columns
%             end
%             tableCondListUnique = unique(tableCondListFull,'stable'); % unique cond names
%             tableStatList{iSheet} = unique(tableStatListFull,'stable'); % unique stat names
            % stat names
            tableStatListFull = dataTables{iSheet}(1,find(~cellfun(@isempty,regexp(dataTables{iSheet}(1,:),regexptranslate('wildcard','*set*'))))+1:end); % list of all stat names in order of columns
            tableStatList{iSheet} = unique(tableStatListFull,'stable'); % unique stat names
            % create column index for each condition on current sheet
%             for iC = 1:length(tableCondListUnique);
%                 condIdx{iC} = find(strcmp(tableCondListFull,tableCondListUnique{iC})) + condStatStrIndex(1) - 1; % create an index for columns for each condition
%             end
            % create column index for each stat on current sheet
            for iS = 1:length(tableStatList{iSheet});
                statIdx{iS} = find(~cellfun(@isempty,regexp(dataTables{iSheet}(1,:),regexptranslate('wildcard',tableStatList{iSheet}{iS})))); % create an index for columns for each stat
            end
            % add conds and indices for current table to condsAll and condsAllIdx lists
%             tableCondList = cellfun(@(c) [sheetsCurrent{iSheet} ':' c],tableCondListUnique,'UniformOutput',false); % label each condition with its spreadsheet name
%             if all(cellfun(@isempty,condsAll))
%                 condsAll = tableCondList';
%             else
%                 condsAll = [condsAll; tableCondList'];
%             end
%             if all(cellfun(@isempty,condsAllIdx))
%                 condsAllIdx = condIdx';
%             else
%                 condsAllIdx = [condsAllIdx; condIdx'];
%             end
%             if isempty(condsAllWhichSheetIdx)
%                 condsAllWhichSheetIdx = repmat(iSheet,[length(tableCondList) 1]);
%             else
%                 condsAllWhichSheetIdx = [condsAllWhichSheetIdx; repmat(iSheet,[length(tableCondList) 1])];
%             end
            % add stat indices for current table to statsAllIdx list
            statsAllIdx{iSheet} = statIdx;
        end
        % store subject codes and stat names subjectsAll and statsAll list
        if iSheet > 1
            if isequal(tableSubjects{1:iSheet}) % use first list if all are identical
                subjectsAll = tableSubjects{1};
            else % if not, use the shorter list of subjects
                warning('The selected sheets do not contain all the same subjects.')
                uselistIdx = cellfun(@length,tableSubjects)==min(cellfun(@length,tableSubjects));
                subjectsAll = tableSubjects{uselistIdx};
                % check that there are no subjects in list that are not in all sheets
                checklistIdx = find(uselistIdx==0);
                for iCheck = 1:length(checklistIdx)
                    subjectsAll = intersect(subjectsAll,tableSubjects{checklistIdx(iCheck)});
                end
            end
            if isequal(tableStatList{1:iSheet})
                statsAll = tableStatList{1};
            else
                warning('The selected sheets do not contain all the same statistics.')
            end
        else % if there is only one current sheet
            subjectsAll = tableSubjects{1};
            statsAll = tableStatList{1};
        end
        % use all current subjects and condtions unless otherwise designated
        subjectsCurrent = subjectsAll;
%         condsCurrent = condsAll;
    end

    function Spreadsheets_CB(~,~)
        % select spreadsheets
        sheetsAll = {'MT Loc','Contrast','Suppression','Summation','V1 Loc','V1_fix Loc','ftap'};
        if all(cellfun(@isempty,sheetsCurrent))
            sheetsSelected = sheetsAll(listdlg('liststring',sheetsAll,'PromptString','Select Experiment(s)'));
            sheetsCurrent = sheetsSelected;
            % reset other variables to prior values
%             condsCurrent = {''}; % to pass between callback functions
            statsCurrent = {''};
            subjectsCurrent = {''};
            filterIdx = {[]};
            filtstatSelected = '';
            operatorSelected = '';
            criterionSelected = [];
            highlight_subjects = {''};
            dataTables = {[]}; % to be defined in LoadData
%             isBlockSheet = [];
%             condsAll = {''};
%             condsAllIdx = {[]};
%             condsAllWhichSheetIdx = [];
            statsAll = {''};
            statsAllIdx = {[]};
            subjectsAll = {''};
            if ~all(cellfun(@isempty,sheetsCurrent))
                LoadData();
            end
        else
            sheetsSelected = sheetsAll(listdlg('liststring',sheetsAll,'initialvalue',find(ismember(sheetsAll,sheetsCurrent)),'PromptString','Select Spreadsheet(s)'));
            if numel(sheetsCurrent) ~= numel(sheetsSelected) || ~all(ismember(sheetsCurrent,sheetsSelected)) || ~all(ismember(sheetsSelected,sheetsCurrent))
                sheetsCurrent = sheetsSelected;
                % reset other variables to prior values
%                 condsCurrent = {''}; % to pass between callback functions
                statsCurrent = {''};
                subjectsCurrent = {''};
                filterIdx = {''};
                filtstatSelected = '';
                operatorSelected = '';
                criterionSelected = [];
                highlight_subjects = {''};
                dataTables = {[]}; % to be defined in LoadData
%                 isBlockSheet = [];
%                 condsAll = {''};
%                 condsAllIdx = {[]};
%                 condsAllWhichSheetIdx = [];
                statsAll = {''};
                statsAllIdx = {[]};
                subjectsAll = {''};
                if ~all(cellfun(@isempty,sheetsCurrent))
                    LoadData();
                end
            end
        end
    end

%     function Conditions_CB(~,~)
%         % must run spreadsheets if condsAll is empty
%         if all(cellfun(@isempty,condsAll))
%             Spreadsheets_CB();
%         end
%         % select conditions
%         if all(cellfun(@isempty,condsCurrent))
%             condsSelected = condsAll(listdlg('liststring',condsAll,'PromptString','Select Condition(s)'));
%             condsCurrent = condsSelected; 
%             if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
%                 RePlot(); 
%             end
%         else
%             condsSelected = condsAll(listdlg('liststring',condsAll,'initialvalue',find(ismember(condsAll,condsCurrent)),'PromptString','Select Condition(s)'));
%             % replot?
%             if numel(condsCurrent) ~= numel(condsSelected) || ~all(ismember(condsCurrent,condsSelected)) || ~all(ismember(condsSelected,condsCurrent))
%                condsCurrent = condsSelected; 
%                if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
%                     RePlot(); 
%                 end 
%             end
%         end
%     end

    function Statistics_CB(~,~)
         % must run spreadsheets if statsAll is empty
        if all(cellfun(@isempty,statsAll))
            Spreadsheets_CB();
        end
        % select statistics
        if all(cellfun(@isempty,statsCurrent))
            statsSelected = statsAll(listdlg('liststring',statsAll,'PromptString','Select Statistic(s)'));
            statsCurrent = statsSelected;
%             if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
            if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                RePlot(); 
            end 
        else
            statsSelected = statsAll(listdlg('liststring',statsAll,'initialvalue',find(ismember(statsAll,statsCurrent)),'PromptString','Select Statistic(s)'));
            % replot?
            if numel(statsCurrent) ~= numel(statsSelected) || ~all(ismember(statsCurrent,statsSelected)) || ~all(ismember(statsSelected,statsCurrent))
               statsCurrent = statsSelected; 
%                if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                    RePlot(); 
               end 
            end
        end
    end
    
    function Subjects_CB(~,~)
        % must run spreadsheets if subjectsAll is empty
        if all(cellfun(@isempty,subjectsAll))
            Spreadsheets_CB();
        end
        % select statistics
        if all(cellfun(@isempty,subjectsCurrent))
            subjectsSelected = subjectsAll(listdlg('liststring',subjectsAll,'PromptString','Select Subject(s)'));
            subjectsCurrent = subjectsSelected;
            highlight_subjects = {[]}; % reset highlighted subjects
%             if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
            if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                RePlot(); 
            end 
        else
            subjectsSelected = subjectsAll(listdlg('liststring',subjectsAll,'initialvalue',find(ismember(subjectsAll,subjectsCurrent)),'PromptString','Select Subject(s)'));
            % replot?
            if numel(subjectsCurrent) ~= numel(subjectsSelected) || ~all(ismember(subjectsCurrent,subjectsSelected)) || ~all(ismember(subjectsSelected,subjectsCurrent))
               subjectsCurrent = subjectsSelected;
               highlight_subjects = {[]}; % reset highlighted subjects
%                if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
               if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                    RePlot(); 
               end 
            end
        end
    end

    function Filter_CB(~,~)
        % select filter type
        filterOptions = {'Clear Filters','Filter by Subject(s)','Filter by Stat(s)'};
        filterSelected = filterOptions(listdlg('liststring',filterOptions,'SelectionMode','single','PromptString','Select from Filtering Options'));
        % designate filter
        if strcmp(filterSelected,filterOptions{1}); % remove filters
            filterIdx = {[]};
            highlight_subjects = {''};
%             if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
            if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                RePlot(); 
            end 
        elseif strcmp(filterSelected,filterOptions{2}); % go to HighlightSubjects for filtering by subject
            HighlightSubjects();
        elseif strcmp(filterSelected,filterOptions{3}); % go to StatFilter for filtering by stats
            StatFilter();
        end
    end
    
    function HighlightSubjects()
        % must run spreadsheets if subjectsAll is empty
        if all(cellfun(@isempty,subjectsAll))
            Spreadsheets_CB();
        end
        % must run subjects if subjectsCurrent is empty
        if all(cellfun(@isempty,subjectsCurrent))
            Subjects_CB();
        end
        % select highlighted subjects
        if all(cellfun(@isempty,highlight_subjects))
            subjectsSelected = subjectsCurrent(listdlg('liststring',subjectsCurrent,'PromptString','Select Subject(s) to Highlight'));
            highlight_subjects = subjectsSelected;
%             if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
            if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                RePlot(); 
            end 
        else
            subjectsSelected = subjectsCurrent(listdlg('liststring',subjectsCurrent,'initialvalue',find(ismember(subjectsCurrent,highlight_subjects)),'PromptString','Select Subject(s) to Highlight'));
            % replot?
            if numel(highlight_subjects) ~= numel(subjectsSelected) || ~all(ismember(highlight_subjects,subjectsSelected)) || ~all(ismember(subjectsSelected,highlight_subjects))
                highlight_subjects = subjectsSelected; 
%                 if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent)) 
                    RePlot(); 
                end 
            end
        end
    end

    function StatFilter()
        % must run spreadsheets if statsAll is empty
        if all(cellfun(@isempty,statsAll))
            Spreadsheets_CB();
        end
        
        % select filter based on stat, logical operator and criterion:
        % create figure
        filtH = figure('tag','Designate Filter by Stat', 'name','Designate Filter by Stat','toolbar','none','numbertitle','off');
        % get figure position and size
        figPosition = getpixelposition(filtH);
        xRange = figPosition(3);
        yRange = figPosition(4);
        % select stat
        if ~exist('statSelectedValues','var')
            statSelectedValues = {[]}; % start empty
        end
        if ~isempty(filtstatSelected)
            uicontrol('style','popup','string',statsAll,'value',find(strcmp(statsAll,filtstatSelected)),'callback',@AssignFilterStat_CB,'pos',[xRange/2-75 yRange*3/4 150 20]);
        else
            uicontrol('style','popup','string',statsAll,'callback',@AssignFilterStat_CB,'pos',[xRange/2-75 yRange*3/4 150 20]);
        end
        uicontrol('style','text','string','Select Statistic by which to Filter','pos',[xRange/2-100 yRange*3/4+30 200 20]);
        % select logical operator
        if ~exist('operatorSelected','var')
            operatorSelected = ''; % start empty
        end
        logicalOperators = {'==','>','<','>=','<='};
        if ~isempty(operatorSelected)
            uicontrol('style','popup','string',logicalOperators,'value',find(strcmp(logicalOperators,operatorSelected)),'callback',@AssignFilterOperator_CB,'pos',[xRange/2-40 yRange*2/4 80 20]);
        else
            uicontrol('style','popup','string',logicalOperators,'callback',@AssignFilterOperator_CB,'pos',[xRange/2-40 yRange*2/4 80 20]);
        end
        uicontrol('style','text','string','Select Logical Operator','pos',[xRange/2-100 yRange*2/4+30 200 20]);
        % enter numarical criterion
        if ~exist('criterionSelected','var')
            criterionSelected = []; % start empty
        end
        if ~isempty(criterionSelected)
            uicontrol('style','edit','string',num2str(criterionSelected),'callback',@AssignFilterCriterion_CB,'pos',[xRange/2-75 yRange*1/4 150 20]);
        else
            uicontrol('style','edit','callback',@AssignFilterCriterion_CB,'pos',[xRange/2-75 yRange*1/4 150 20]);
        end
        uicontrol('style','text','string','Select Numerical Criterion','pos',[xRange/2-100 yRange*1/4+30 200 20]);
        % apply filter:
        uicontrol('style','pushbutton','string','Apply Filter','callback',@ApplyFilter_CB,'pos',[xRange/2-40 50 80 30]);
        
        % callback functions for filtering
        function AssignFilterStat_CB(source,~)
            clear filtstatSelected statSelectedValues
            % assign statSelected
            index = source.Value;
            list = source.String;
            filtstatSelected = list{index};
            % get values of stat for each sheet
            for iSheet = 1:length(sheetsCurrent)
                clear statColumnIdx
                statColumnIdx = ~cellfun(@isempty,regexp(dataTables{iSheet}(1,:),regexptranslate('wildcard',filtstatSelected)));
                statSelectedValues{iSheet} = dataTables{iSheet}(:,statColumnIdx);
            end
        end
        function AssignFilterOperator_CB(source,~)
            clear operatorSelected
            % assign statSelected
            index = source.Value;
            list = source.String;
            operatorSelected = list{index};
        end
        function AssignFilterCriterion_CB(source,~)
            clear criterionSelected
            % assign statSelected
            criterionSelected = str2double(source.String);
        end
        function ApplyFilter_CB(~,~)
            if exist('statSelectedValues','var') && ~all(cellfun(@isempty,statSelectedValues)) && exist('operatorSelected','var') && ischar(operatorSelected) && ~isempty(operatorSelected) && exist('criterionSelected','var') && ~isnan(criterionSelected) && ~isempty(criterionSelected)
                % generate filterIdx per sheet in sheetsCurrent and per cond in each sheet
                for iS = 1:length(statSelectedValues)
                    clear thisSheetStats
                    thisSheetStats = statSelectedValues{iS}; % extract data
                    thisSheetStats(~cellfun(@isnumeric,thisSheetStats)) = {nan}; % replace column labels
                    thisSheetStats = cell2mat(thisSheetStats); % convert to matrix
                    filterIdx(iS) = {eval([mat2str(thisSheetStats) operatorSelected num2str(criterionSelected)])}; % evaluate filter for all values of matrix
                end
                % close filter selection interface
                close gcf
                clear filtH
                % replot
%                 if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                    RePlot(); 
                end
            else
                warning('Must select filter parameters before applying filter');
                % start over with filter creation
                close gcf
                clear filtH
                StatFilter();
            end
        end    
    end

    function RePlot()
        % preset counter
        x_labelIdx = 1;
        % separate into vectors of data by group, stats of interest, generate labels, and group information for plotting
        for iSheet = 1:length(sheetsCurrent);
            % indices for group membership
            ASDindex{iSheet} = find(~cellfun(@isempty,regexp(dataTables{iSheet}(:,1),regexptranslate('wildcard','*G1*')))==1);
            NTindex{iSheet} = find(~cellfun(@isempty,regexp(dataTables{iSheet}(:,1),regexptranslate('wildcard','*G3*')))==1);
            % use only subjects currently selected
            if ~all(cellfun(@isempty,subjectsCurrent))
                ASDindex{iSheet} = ASDindex{iSheet}(ismember(dataTables{iSheet}(ASDindex{iSheet},1),subjectsCurrent));
                NTindex{iSheet} = NTindex{iSheet}(ismember(dataTables{iSheet}(NTindex{iSheet},1),subjectsCurrent));
            end
            % index for subjects of interest (highlight_subjects)
            if ~all(cellfun(@isempty,highlight_subjects))
                ASDhighlightIndex{iSheet} = zeros(length(ASDindex{iSheet}),1);
                NThighlightIndex{iSheet} = zeros(length(NTindex{iSheet}),1);
                for iSubj = 1:length(highlight_subjects);
                    ASDhighlightIndex{iSheet} = ASDhighlightIndex{iSheet}+strcmp(dataTables{iSheet}(ASDindex{iSheet},1),highlight_subjects{iSubj});
                    NThighlightIndex{iSheet} = NThighlightIndex{iSheet}+strcmp(dataTables{iSheet}(NTindex{iSheet},1),highlight_subjects{iSubj});
                end
                ASDhighlightIndex{iSheet} = logical(ASDhighlightIndex{iSheet});
                NThighlightIndex{iSheet} = logical(NThighlightIndex{iSheet});
            end
%         end
        
        %%% MAJOR CHANGE IN INDEXING (2/7/17)
        % in the following systax, indexing by condition has been replaced
        % by indexing by sheet, where each sheet represents an experiment
        
        % assignment of data arrays
%         for iCond = 1:length(sheetsCurrent)
%             % index for this cond and the associated sheet
%             clear thisCondIdx thisCondSheetIdx
%             thisCondIdx = strcmp(condsAll,condsCurrent{iCond}); % idx in condsAll associated with iCond
%             thisCondSheetIdx = condsAllWhichSheetIdx(strcmp(condsAll,condsCurrent{iCond})); % sheet idx associated with iCond
            for iStat = 1:length(statsCurrent)
                % index for this stat 
                clear otherStats thisStatIdx
                otherStats = setdiff(statsAll,statsCurrent{iStat});
                thisStatIdx = strcmp(statsAll,statsCurrent{iStat}) & ~cellfun(@(x) any(strcmp(otherStats,x)),statsAll); % idx in statsAll associated with iStat
                % store column values by stat and sheet
                clear statSheetIdx
                statSheetIdx = statsAllIdx{iSheet}{thisStatIdx}; % should be only one column per stat per sheet
                % assign data to plot to stat structure
                stat(iStat).name = statsCurrent{iStat};
                stat(iStat).ASDvalues{iSheet} = arrayfun(@cell2mat,dataTables{iSheet}(ASDindex{iSheet},statSheetIdx)); % for ASD group
                stat(iStat).NTvalues{iSheet} = arrayfun(@cell2mat,dataTables{iSheet}(NTindex{iSheet},statSheetIdx)); % for NT group
            
                % create mouseover labels for subject,fMRI session, and
                % set# and prepare exportData matrix
%                 switch isBlockSheet(thisCondSheetIdx) % labels different for block and summary stats
%                     case 0 % for summary stats:
                for iA = 1:length(ASDindex{iSheet})
                    stat(iStat).ASDlabels{iSheet}{iA} = [dataTables{iSheet}{ASDindex{iSheet}(iA),1} ': ' dataTables{iSheet}{ASDindex{iSheet}(iA),2} ': set' num2str(dataTables{iSheet}{ASDindex{iSheet}(iA),3})];
                end
                for iT = 1:length(NTindex{iSheet})
                    stat(iStat).NTlabels{iSheet}{iT} = [dataTables{iSheet}{NTindex{iSheet}(iT),1} ': ' dataTables{iSheet}{NTindex{iSheet}(iT),2} ': set' num2str(dataTables{iSheet}{NTindex{iSheet}(iT),3})];
                end

                %%% format data for export (01/26/17):
                for iSubj = 1:length(subjectsCurrent)
                    for iSess = 1:2
                        for iSet = 1:9
                            % indices for this subject, session, and set
                            if isfield(stat,'ASDlabels') && isfield(stat,'NTlabels')
                                thisSubjIdx = ~cellfun(@isempty,regexp([stat(iStat).ASDlabels{iSheet} stat(iStat).NTlabels{iSheet}],regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
                                thisSessIdx = ~cellfun(@isempty,regexp([stat(iStat).ASDlabels{iSheet} stat(iStat).NTlabels{iSheet}],regexptranslate('wildcard',['*fMRI' num2str(iSess) '*'])));
                                thisSetIdx = ~cellfun(@isempty,regexp([stat(iStat).ASDlabels{iSheet} stat(iStat).NTlabels{iSheet}],regexptranslate('wildcard',['*set' num2str(iSet) '*'])));
                            elseif isfield(stat,'ASDlabels') && ~isfield(stat,'NTlabels')
                                thisSubjIdx = ~cellfun(@isempty,regexp(stat(iStat).ASDlabels{iSheet},regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
                                thisSessIdx = ~cellfun(@isempty,regexp(stat(iStat).ASDlabels{iSheet},regexptranslate('wildcard',['*fMRI' num2str(iSess) '*'])));
                                thisSetIdx = ~cellfun(@isempty,regexp(stat(iStat).ASDlabels{iSheet},regexptranslate('wildcard',['*set' num2str(iSet) '*'])));
                            elseif ~isfield(stat,'ASDlabels') && isfield(stat,'NTlabels')
                                thisSubjIdx = ~cellfun(@isempty,regexp(stat(iStat).NTlabels{iSheet},regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
                                thisSessIdx = ~cellfun(@isempty,regexp(stat(iStat).NTlabels{iSheet},regexptranslate('wildcard',['*fMRI' num2str(iSess) '*'])));
                                thisSetIdx = ~cellfun(@isempty,regexp(stat(iStat).NTlabels{iSheet},regexptranslate('wildcard',['*set' num2str(iSet) '*'])));
                            elseif ~isfield(stat,'ASDlabels') && ~isfield(stat,'NTlabels')
                                warning('No subjects of either group selected');
                                thisSubjIdx = false;
                                thisSessIdx = false;
                                thisSetIdx = false;
                            end
                            % store data
                            tempData = [stat(iStat).ASDvalues{iSheet}; stat(iStat).NTvalues{iSheet}];
                            if isempty(find(thisSubjIdx & thisSessIdx & thisSetIdx,1))
                                exportData(iSubj,iSess,iSet,iSheet,iStat) = nan;
                            else
                                exportData(iSubj,iSess,iSet,iSheet,iStat) = tempData(thisSubjIdx & thisSessIdx & thisSetIdx);
                            end
                        end
                    end
                end
                        
%                     case 1 % for block stats:
%                         maxBlocks = 0; 
%                         for iA = 1:length(ASDindex{thisCondSheetIdx})
%                             % for plots
%                             stat(iStat).ASDlabels{iCond}{iA} = [dataTables{thisCondSheetIdx}{ASDindex{thisCondSheetIdx}(iA),1} ': ' dataTables{thisCondSheetIdx}{ASDindex{thisCondSheetIdx}(iA),2} ': ' dataTables{thisCondSheetIdx}{ASDindex{thisCondSheetIdx}(iA),3} ': block' num2str(dataTables{thisCondSheetIdx}{ASDindex{thisCondSheetIdx}(iA),5})];
%                             % for export
%                             if maxBlocks < dataTables{thisCondSheetIdx}{ASDindex{thisCondSheetIdx}(iA),5}
%                                 maxBlocks = dataTables{thisCondSheetIdx}{ASDindex{thisCondSheetIdx}(iA),5};
%                             end
%                         end
%                         for iT = 1:length(NTindex{thisCondSheetIdx})
%                             % for plots
%                             stat(iStat).NTlabels{iCond}{iT} = [dataTables{thisCondSheetIdx}{NTindex{thisCondSheetIdx}(iT),1} ': ' dataTables{thisCondSheetIdx}{NTindex{thisCondSheetIdx}(iT),2} ': ' dataTables{thisCondSheetIdx}{NTindex{thisCondSheetIdx}(iT),3} ': block' num2str(dataTables{thisCondSheetIdx}{NTindex{thisCondSheetIdx}(iT),5})];
%                             % for export
%                             if maxBlocks < dataTables{thisCondSheetIdx}{NTindex{thisCondSheetIdx}(iT),5}
%                                 maxBlocks = dataTables{thisCondSheetIdx}{NTindex{thisCondSheetIdx}(iT),5};
%                             end
%                         end
%                         
%                         %%% format data for export (01/26/17):
%                         for iSubj = 1:length(subjectsCurrent)
%                             for iSess = 1:2
%                                 for iSet = 1:9
%                                     for iBlock = 1:maxBlocks
%                                         % indices for this subject, session, and set
%                                         thisSubjIdx = ~cellfun(@isempty,regexp([stat(iStat).ASDlabels{iCond} stat(iStat).NTlabels{iCond}],regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
%                                         thisSessIdx = ~cellfun(@isempty,regexp([stat(iStat).ASDlabels{iCond} stat(iStat).NTlabels{iCond}],regexptranslate('wildcard',['*fMRI' num2str(iSess) '*'])));
%                                         thisSetIdx = ~cellfun(@isempty,regexp([stat(iStat).ASDlabels{iCond} stat(iStat).NTlabels{iCond}],regexptranslate('wildcard',['*set' num2str(iSet) '*'])));
%                                         thisBlockIdx = ~cellfun(@isempty,regexp([stat(iStat).ASDlabels{iCond} stat(iStat).NTlabels{iCond}],regexptranslate('wildcard',['*block' num2str(iBlock) '*'])));
%                                         % store data
%                                         tempData = [stat(iStat).ASDvalues{iCond}; stat(iStat).NTvalues{iCond}];
%                                         if isempty(find(thisSubjIdx & thisSessIdx & thisSetIdx & thisBlockIdx,1))
%                                             exportData(iSubj,iSess,iSet,iCond,iStat) = nan;
%                                         else
%                                             exportData(iSubj,iSess,iSet,iBlock,iCond,iStat) = tempData(thisSubjIdx & thisSessIdx & thisSetIdx & thisBlockIdx);
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                 end
            end
            
            % use filters only for selected conds 
%             if ~all(cellfun(@isempty,filterIdx))
%                 % generate index for the condition number within the current sheet for accessing filterIdx (easier to do this here rather than later)
%                 clear thisSheetCondIdx
%                 thisSheetCondIdx = strcmp(condsAll(~cellfun(@isempty,regexp(condsAll,regexptranslate('wildcard',[sheetsCurrent{thisCondSheetIdx} ':'])))),condsCurrent{iCond});
%                 % set up useFiltIdx as a cell array with one cell per iCond
%                 % this will allow for easier indexing in the creation of highlightIndices
%                 useFiltIdx{iCond} = filterIdx{thisCondSheetIdx}(:,thisSheetCondIdx);
%             end
%         end

        % x-axis labels
%         x_labelIdx = 1;
%         for iCond = 1:length(condsCurrent)
            % do this differently depending on how selected subjects are grouped
            if any(cellfun(@any,ASDindex)) && any(cellfun(@any,NTindex)) % if there are some subjects of each group
                x_labels{x_labelIdx} = [sheetsCurrent{iSheet} ' ASD']; % ASD label
                x_labels{x_labelIdx+1} = [sheetsCurrent{iSheet} ' NT']; % NT label
                x_labels{x_labelIdx+2} = ''; % blank label for spacing
                x_labelIdx = length(x_labels)+1; % advance counter
            elseif any(cellfun(@any,ASDindex)) && ~any(cellfun(@any,NTindex)) % if there are only ASD subjects
                x_labels{x_labelIdx} = [sheetsCurrent{iSheet} ' ASD']; % ASD label
                x_labels{x_labelIdx+1} = ''; % blank label for spacing
                x_labelIdx = length(x_labels)+1; % advance counter
            elseif ~any(cellfun(@any,ASDindex)) && any(cellfun(@any,NTindex)) % if there are only NT subjects
                x_labels{x_labelIdx} = [sheetsCurrent{iSheet} ' NT']; % ASD label
                x_labels{x_labelIdx+1} = ''; % blank label for spacing
                x_labelIdx = length(x_labels)+1; % advance counter
            elseif ~any(cellfun(@any,ASDindex)) && ~any(cellfun(@any,NTindex)) % if there are no subjects of either group
                x_labels{x_labelIdx} = ''; % blank label for spacing
                x_labelIdx = length(x_labels)+1; % advance counter
            end
        end
        x_labels(end) = []; % remove last blank x tick
        
        
        % group data arrays, mouseover labels, and highlight indices for plotting
        if ~all(cellfun(@isempty,highlight_subjects)) || ~all(cellfun(@isempty,filterIdx))
            highlightIndices = [];
            distIdx = []; % index into highlightIndices by x_label number
        end
        % setup x axis label
        xAxisLabel = {}; % to contain one name per experiment/spreadsheet
        xAxisLabelX = []; % to contain x coordinates at which to place experiment labels
        for iStat = 1:length(statsCurrent)
           for iX = 1:length(x_labels)
               if isempty(x_labels{iX})
                   % match empty x_labes with empty data arrays and empty datapoint labels 
                   stat(iStat).plotData(iX) = {nan(1)};
                   stat(iStat).plotLabels(iX) = {''};
                   stat(iStat).x_labels(iX) = x_labels(iX);
                   % append highlightIndices
                   if ( ~all(cellfun(@isempty,highlight_subjects)) || ~all(cellfun(@isempty,filterIdx)) ) && iStat==1
                       highlightIndices = vertcat(highlightIndices,0);
                       distIdx = vertcat(distIdx,iX);
                   end
               else 
                   % indices for condsCurrent and sheetsCurrent for this x_label:
                   clear condLabelIdx sheetLabelIdx
                   % do this differently for block and non-block x_labels
%                    if ~isempty(strfind(x_labels{iX},'block')) % if this x_label is for stats by block
%                        condLabelIdx = find(AK_whichPattern(x_labels{iX},condsCurrent) & ~cellfun(@isempty,strfind(condsCurrent,'block'))'); % which cond corresponds to x_labels{iX}
%                        sheetLabelIdx = find(AK_whichPattern(x_labels{iX},sheetsCurrent) & ~cellfun(@isempty,strfind(sheetsCurrent,'block'))); % which sheet corresponds to x_labels{iX}
%                    else % if this x_label is for stats by set, rather than by block
%                        condLabelIdx = find(AK_whichPattern(x_labels{iX},condsCurrent) & cellfun(@isempty,strfind(condsCurrent,'block'))'); % which cond corresponds to x_labels{iX}
%                        sheetLabelIdx = find(AK_whichPattern(x_labels{iX},sheetsCurrent) & cellfun(@isempty,strfind(sheetsCurrent,'block'))); % which sheet corresponds to x_labels{iX}
%                    end
%                    condLabelIdx = find(AK_whichPattern(x_labels{iX},condsCurrent)); % which cond corresponds to x_labels{iX}
                   sheetLabelIdx = find(AK_whichPattern(x_labels{iX},sheetsCurrent)); % which sheet corresponds to x_labels{iX} 
                   if find(AK_whichPattern(x_labels{iX},{'ASD','NT'})) == 1 && ~cellfun(@isempty,stat(iStat).ASDvalues(sheetLabelIdx)) % is ASD label?
                       % match ASD x_labels for cond with proper data arrays and empty datapoint labels
                       stat(iStat).plotData(iX) = stat(iStat).ASDvalues(sheetLabelIdx);
                       stat(iStat).plotLabels(iX) = stat(iStat).ASDlabels(sheetLabelIdx);
                       stat(iStat).x_labels(iX) = {'ASD'};
                       % append highlightIndices
                       if ~all(cellfun(@isempty,highlight_subjects)) && all(cellfun(@isempty,filterIdx)) && iStat==1
                           highlightIndices = vertcat(highlightIndices,ASDhighlightIndex{sheetLabelIdx});
                           distIdx = vertcat(distIdx,repmat(iX,size(ASDhighlightIndex{sheetLabelIdx})));
                       end
                       if all(cellfun(@isempty,highlight_subjects)) && ~all(cellfun(@isempty,filterIdx)) && iStat==1
%                            highlightIndices = vertcat(highlightIndices,2.*useFiltIdx{condLabelIdx}(ASDindex{sheetLabelIdx}));
%                            distIdx = vertcat(distIdx,repmat(iX,size(2.*useFiltIdx{condLabelIdx}(ASDindex{sheetLabelIdx}))));
                           highlightIndices = vertcat(highlightIndices,2.*filterIdx{sheetLabelIdx}(ASDindex{sheetLabelIdx}));
                           distIdx = vertcat(distIdx,repmat(iX,size(2.*filterIdx{sheetLabelIdx}(ASDindex{sheetLabelIdx})))); 
                       end
                       if ~all(cellfun(@isempty,highlight_subjects)) && ~all(cellfun(@isempty,filterIdx)) && iStat==1
%                            highlightIndices = vertcat(highlightIndices, ASDhighlightIndex{sheetLabelIdx} + 2.*useFiltIdx{condLabelIdx}(ASDindex{sheetLabelIdx}) );
%                            distIdx = vertcat(distIdx,repmat(iX,size(ASDhighlightIndex{sheetLabelIdx} + 2.*useFiltIdx{condLabelIdx}(ASDindex{sheetLabelIdx}))));
                           highlightIndices = vertcat(highlightIndices, ASDhighlightIndex{sheetLabelIdx} + 2.*filterIdx{sheetLabelIdx}(ASDindex{sheetLabelIdx}) );
                           distIdx = vertcat(distIdx,repmat(iX,size(ASDhighlightIndex{sheetLabelIdx} + 2.*filterIdx{sheetLabelIdx}(ASDindex{sheetLabelIdx}))));
                       end
                       % append x axis label
                       if iStat==1
                           xAxisLabel = [xAxisLabel {x_labels{iX}(1:end-4)}]; % add experiment name to x axis label
                           xAxisLabelX = [xAxisLabelX iX];
                       end
                   elseif find(AK_whichPattern(x_labels{iX},{'ASD','NT'})) == 2 && ~cellfun(@isempty,stat(iStat).NTvalues(sheetLabelIdx)) % is NT label?
                       % match NT x_labels for cond with proper data arrays and empty datapoint labels
                       stat(iStat).plotData(iX) = stat(iStat).NTvalues(sheetLabelIdx);
                       stat(iStat).plotLabels(iX) = stat(iStat).NTlabels(sheetLabelIdx);
                       stat(iStat).x_labels(iX) = {'NT'};
                       % append highlightIndices
                       if ~all(cellfun(@isempty,highlight_subjects)) && all(cellfun(@isempty,filterIdx)) && iStat==1
                           highlightIndices = vertcat(highlightIndices,NThighlightIndex{sheetLabelIdx});
                           distIdx = vertcat(distIdx,repmat(iX,size(NThighlightIndex{sheetLabelIdx})));
                       end
                       if all(cellfun(@isempty,highlight_subjects)) && ~all(cellfun(@isempty,filterIdx)) && iStat==1
%                            highlightIndices = vertcat(highlightIndices,2.*useFiltIdx{condLabelIdx}(NTindex{sheetLabelIdx}));
%                            distIdx = vertcat(distIdx,repmat(iX,size(2.*useFiltIdx{condLabelIdx}(NTindex{sheetLabelIdx}))));
                           highlightIndices = vertcat(highlightIndices,2.*filterIdx{sheetLabelIdx}(NTindex{sheetLabelIdx}));
                           distIdx = vertcat(distIdx,repmat(iX,size(2.*filterIdx{sheetLabelIdx}(NTindex{sheetLabelIdx}))));
                       end
                       if ~all(cellfun(@isempty,highlight_subjects)) && ~all(cellfun(@isempty,filterIdx)) && iStat==1
%                            highlightIndices = vertcat(highlightIndices, NThighlightIndex{sheetLabelIdx} + 2.*useFiltIdx{condLabelIdx}(NTindex{sheetLabelIdx}) );
%                            distIdx = vertcat(distIdx,repmat(iX,size(NThighlightIndex{sheetLabelIdx} + 2.*useFiltIdx{condLabelIdx}(NTindex{sheetLabelIdx}))));
                           highlightIndices = vertcat(highlightIndices, NThighlightIndex{sheetLabelIdx} + 2.*filterIdx{sheetLabelIdx}(NTindex{sheetLabelIdx}) );
                           distIdx = vertcat(distIdx,repmat(iX,size(NThighlightIndex{sheetLabelIdx} + 2.*filterIdx{sheetLabelIdx}(NTindex{sheetLabelIdx}))));
                       end
                       % append x axis label (if this is the first stat and 
                       % (either there are no ASD x labels or the preceding x label was not an ASD label)) 
                       if iStat==1 && ( ~any(cell2mat(cellfun(@(x) find(AK_whichPattern(x,{'ASD','NT'}))==1,x_labels(~cellfun(@isempty,x_labels)),'UniformOutput',false))) || find(AK_whichPattern(x_labels{iX-1},{'ASD','NT'})) ~= 1 ) % only append label if this sheet name has not already been added to the x axis label
                           xAxisLabel = [xAxisLabel {x_labels{iX}(1:end-3)}]; % add experiment name to x axis label
                           xAxisLabelX = [xAxisLabelX iX];
                       end
                   end
               end
           end
        end

        % figures
        close gcf
        clear fH
        InitializeFigureWindow();
        clear structure
        if ~all(cellfun(@isempty,highlight_subjects)) || ~all(cellfun(@isempty,filterIdx)) % if a filter has been applied
            for iStat = 1:length(statsCurrent)
                % plot
                subplot(length(statsCurrent),1,iStat);
%                 stat(iStat).handles = AK_plotSpread_mouseover(stat(iStat).plotData,stat(iStat).plotLabels,'categoryIdx',highlightIndices,'xNames',stat(iStat).x_labels,'showMM',5);
                stat(iStat).handles = AK_plotSpread_forGUI(stat(iStat).plotData,'categoryIdx',highlightIndices,'xNames',stat(iStat).x_labels,'showMM',5);                
                title(stat(iStat).name);
                ylabel('mm'); % all framewise displacement stats are in units of mm
                set(stat(iStat).handles{2}(2),'LineWidth',1.5)
                % figure out y position of x labels relative to axes
                yL = get(gca,'YLim');
                xAxisLabelY = yL(1) - (0.1 + .05*(length(statsCurrent)-1))*(yL(2)-yL(1)); % lower limit of y axis minus 10% of y axis range, minus an extra 5% per subplot
                % experiment labels as text
                for iE = 1:length(xAxisLabel)
                    text(xAxisLabelX(iE),xAxisLabelY,xAxisLabel{iE});
                end
                % legend
                if iStat==1
                    % prepare legend input:
                    clear subjStr statStr legendCell
                    % subject filter string
                    subjStr = ['Subject(s): ' cellfun(@(x) [x ','],highlight_subjects','UniformOutput',false)]; % cat
                    subjStr{end}(end) = ''; % remove final comma
                    subjStr = strjoin(subjStr); % convert to string
                    % stat filter string
                    statStr = ['Filter: ' filtstatSelected operatorSelected num2str(criterionSelected)];
                    % create legend cell array with one cell per unique highlightIndex
                    legendCell = {'Mean and SD'};
                    hIdxUnique = unique(highlightIndices); % used again later
                    for iU = 1:length(hIdxUnique)
                        if hIdxUnique(iU)==0
                            legendCell = [legendCell,{'Unfiltered'}];
                        elseif hIdxUnique(iU)==1
                            legendCell = [legendCell,{subjStr}];
                        elseif hIdxUnique(iU)==2
                            legendCell = [legendCell,{statStr}];
                        elseif hIdxUnique(iU)==3
                            legendCell = [legendCell,{'Both Filters'}];
                        end
                    end
                    % apply legend
                    legend(gca,legendCell);
                    % OLD LEGEND CODE
%                     if ~all(cellfun(@isempty,highlight_subjects)) && ~all(cellfun(@isempty,filterIdx)) % both kinds of filter
%                         clear subjStr
%                         subjStr = ['Subject(s): ' cellfun(@(x) [x ','],highlight_subjects','UniformOutput',false)]; % cat
%                         subjStr{end}(end) = ''; % remove final comma
%                         subjStr = strjoin(subjStr); % convert to string
%                         if length(unique(highlightIndices))==3 % no overlap in filter coverage
%                             legend(gca,{'Mean and SD','Unfiltered',subjStr,['Filter: ' filtstatSelected operatorSelected num2str(criterionSelected)]});
%                         elseif length(unique(highlightIndices))==4 % some overlap in filter coverage 
%                             legend(gca,{'Mean and SD','Unfiltered',subjStr,['Filter: ' filtstatSelected operatorSelected num2str(criterionSelected)],'Both Filters'});
%                         elseif length(unique(highlightIndices))==2 % no stat met criteria; subjects only
%                             legend(gca,{'Mean and SD','Unfiltered',subjStr});
%                         end
%                     elseif ~all(cellfun(@isempty,highlight_subjects)) && all(cellfun(@isempty,filterIdx)) % subject filter only
%                         clear subjStr
%                         subjStr = ['Subject(s): ' cellfun(@(x) [x ','],highlight_subjects','UniformOutput',false)]; % cat
%                         subjStr{end}(end) = ''; % remove final comma
%                         subjStr = strjoin(subjStr); % convert to string
%                         legend(gca,{'Mean and SD','Unfiltered',subjStr});
%                     elseif all(cellfun(@isempty,highlight_subjects)) && ~all(cellfun(@isempty,filterIdx)) % stat filter only 
%                         legend(gca,{'Mean and SD','Unfiltered',['Filter: ' filtstatSelected operatorSelected num2str(criterionSelected)]});
%                     end
                end
                % prepare output struct (OLD DATA STRUCT)
                plotStruct(iStat).statName = stat(iStat).name;
                plotStruct(iStat).handles = stat(iStat).handles;
                plotStruct(iStat).plotData = stat(iStat).plotData;
                plotStruct(iStat).mouseover_labels = stat(iStat).plotLabels;
                plotStruct(iStat).x_labels = x_labels;
                plotStruct(iStat).filterIdx = highlightIndices; % only where filter has been applied
                plotStruct(iStat).filterKey = legendCell; % only where filter has been applied
            end
        else
            for iStat = 1:length(statsCurrent)
                % plot
                subplot(length(statsCurrent),1,iStat);
%                 stat(iStat).handles = AK_plotSpread_mouseover(stat(iStat).plotData,stat(iStat).plotLabels,'xNames',stat(iStat).x_labels,'showMM',5);
                stat(iStat).handles = AK_plotSpread_forGUI(stat(iStat).plotData,'xNames',stat(iStat).x_labels,'showMM',5);
                title(stat(iStat).name);
                ylabel('mm'); % all framewise displacement stats are in units of mm
                set(stat(iStat).handles{2}(2),'LineWidth',1.5)
                % figure out y position of x labels relative to axes
                yL = get(gca,'YLim');
                xAxisLabelY = yL(1) - (0.1 + .05*(length(statsCurrent)-1))*(yL(2)-yL(1)); % lower limit of y axis minus 10% of y axis range, minus an extra 5% per subplot
                % experiment labels as text
                for iE = 1:length(xAxisLabel)
                    text(xAxisLabelX(iE),xAxisLabelY,xAxisLabel{iE});
                end
                % legend
                if iStat==1
                    legend(gca,{'Mean and SD','Stat Values'});
                end
                % prepare output struct (OLD DATA STRUCT)
                plotStruct(iStat).statName = stat(iStat).name;
                plotStruct(iStat).handles = stat(iStat).handles;
                plotStruct(iStat).plotData = stat(iStat).plotData;
                plotStruct(iStat).mouseover_labels = stat(iStat).plotLabels;
                plotStruct(iStat).x_labels = x_labels;
            end
        end
        % call mouseover function continuously
        if exist('plotStruct','var') && isfield(plotStruct,'plotData') && isfield(plotStruct,'mouseover_labels')
            % prepare inputs to mouseover function as cell arrays with matching dimensions 
            for iSt = 1:length(plotStruct) % stats
                for iDist = 1:length(plotStruct(iSt).plotData) % x_labels
                    % treat filtered data differently b/c of extra dimension added by plotting multiple categories at each x_label
                    if length(plotStruct(iSt).handles{1}(1,:))>1
%                         hIdxUnique = unique(highlightIndices); % assigned when creating legend
                        for iFiltGroup = 1:length(hIdxUnique) % filter groups
                            % assign cells
                            rawData{iSt,iDist,iFiltGroup} = plotStruct(iSt).plotData{iDist}(highlightIndices(distIdx==iDist)==hIdxUnique(iFiltGroup));
                            if ~isempty(plotStruct(iSt).mouseover_labels{iDist})
                                mouseoverLabels{iSt,iDist,iFiltGroup} = plotStruct(iSt).mouseover_labels{iDist}(highlightIndices(distIdx==iDist)==hIdxUnique(iFiltGroup));
                            else
                                mouseoverLabels{iSt,iDist,iFiltGroup} = plotStruct(iSt).mouseover_labels{iDist}; % this syntax avoids indexing an empty array
                            end
                            plotH{iSt,iDist,iFiltGroup} = handle(plotStruct(iSt).handles{1}(iDist,iFiltGroup));
                            % remove nans from both cell arrays at the level of each distribution of points
                            if ~isempty(mouseoverLabels{iSt,iDist,iFiltGroup}) % leave nans marking blank x_labels in place to preserve dimensionality
                                clear badData
                                badData = ~isfinite(rawData{iSt,iDist,iFiltGroup});
                                rawData{iSt,iDist,iFiltGroup}(badData) = [];
                                mouseoverLabels{iSt,iDist,iFiltGroup}(badData) = [];
                            end
                        end
                    else
                        % assign cells
                        rawData{iSt,iDist} = plotStruct(iSt).plotData{iDist};
                        mouseoverLabels{iSt,iDist} = plotStruct(iSt).mouseover_labels{iDist};
                        plotH{iSt,iDist} = handle(plotStruct(iSt).handles{1}(iDist));
                        % remove nans from both cell arrays at the level of each distribution of points
                        if ~isempty(mouseoverLabels{iSt,iDist}) % leave nans marking blank x_labels in place to preserve dimensionality
                            clear badData
                            badData = ~isfinite(rawData{iSt,iDist});
                            rawData{iSt,iDist}(badData) = [];
                            mouseoverLabels{iSt,iDist}(badData) = [];
                        end
                    end
                end
            end
            % pass input to mouseover function
            set(fH,'windowbuttonmotionfcn',{@mouseover,plotH,rawData,mouseoverLabels});
            % pass input to resize function
            set(fH,'ResizeFcn',{@drawXLabels,xAxisLabel,xAxisLabelX})
        end   
    end
    
    function mouseover(src,ev,L,data,labels)
        %since this is a figure callback, the first input is the figure handle:
        f=handle(src);
        %like all callbacks, the second input, ev, isn't used. 

        %determine which object is below the cursor:
        obj=hittest(f); %<-- the important line in this function
        
        if ~isequal(obj,f) % if the current object is not the figure window
            if any(any(any(cellfun(@(x) isequal(obj,x),L)))) %if over any plot... 
                % change pointer when hovering over plot
                set(f,'Pointer','crosshair');

                % determine which stat and distribution plot
                Lindex = find(cellfun(@(x) isequal(obj,x),L));
                [~, Lcol, ~] = ind2sub(size(L),Lindex);

                %get cursor coordinates in its axes:
                a=get(L{Lindex},'parent');
                point=get(a,'currentpoint');
                xclick=point(1,1,1);
                yclick=point(1,2,1);

                %determine which point we're over:
                idx=findclosestpoint2D(xclick,yclick,L{Lindex});

                %make a "tool tip" that displays data point
                xl = xlim; yl = ylim;
                xrange = xl(2)-xl(1); yrange = yl(2)-yl(1);
                xoffset=xrange/1000;
                yoffset=yrange/1000;

                delete(findobj(f,'tag','mytooltip')); %delete last tool tip
                text(Lcol+xoffset,data{Lindex}(idx)+yoffset,labels{Lindex}{idx},...
                    'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
                    'hittest','off');
            else
                % change pointer back and delete plot
                set(gcf,'Pointer','arrow');
                delete(findobj(f,'tag','mytooltip')); %delete last tool tip
            end
        else
            % change pointer back and delete plot
            set(gcf,'Pointer','arrow');
            delete(findobj(f,'tag','mytooltip')); %delete last tool tip
        end
    end

    function index=findclosestpoint2D(xclick,yclick,datasource)
        %this function checks which point in the plotted line "datasource"
        %is closest to the point specified by xclick/yclick. It's kind of 
        %complicated, but this isn't really what this demo is about...

        xdata=get(datasource,'xdata');
        ydata=get(datasource,'ydata');

        activegraph=get(datasource,'parent');

        pos=getpixelposition(activegraph);
        xlim=get(activegraph,'xlim');
        ylim=get(activegraph,'ylim');

        %make conversion factors, units to pixels:
        xconvert=(xlim(2)-xlim(1))/pos(3);
        yconvert=(ylim(2)-ylim(1))/pos(4);

        Xclick=(xclick-xlim(1))/xconvert;
        Yclick=(yclick-ylim(1))/yconvert;

        Xdata=(xdata-xlim(1))/xconvert;
        Ydata=(ydata-ylim(1))/yconvert;

        Xdiff=Xdata-Xclick;
        Ydiff=Ydata-Yclick;

        distnce=sqrt(Xdiff.^2+Ydiff.^2);

        index=distnce==min(distnce);

        index=index(:); %make sure it's a column.

        if sum(index)>1
            thispoint=find(distnce==min(distnce),1);
            index=false(size(distnce));
            index(thispoint)=true;
        end
    end

    function drawXLabels(src,ev,xAxisLabel,xAxisLabelX)
        %since this is a figure callback, the first input is the figure handle:
        f = handle(src);
        a = findobj(f,'Type','Axes');
        % cycle through axes
        for iA = 1:length(a)
            % set current axis
            axes(a(iA));
            % figure out y position of x labels relative to axes
            yL = a(iA).YLim;
            xAxisLabelY = yL(1) - (0.1 + .05*(length(a)-1))*(yL(2)-yL(1)); % lower limit of y axis minus 10% of y axis range, minus an extra 5% per subplot
            % get text handles for current axis
            t = findobj(a(iA),'Type','Text');
            % experiment labels as text
            for iE = 1:length(xAxisLabel)
                set(t(iE),'Visible','off');
                t(iE) = text(xAxisLabelX(iE),xAxisLabelY,xAxisLabel{iE});
                set(t(iE),'Visible','on');
            end
        end
    end

end