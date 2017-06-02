function AK_GABAinASD_EyeTrackingFiguresGUI( xlsfile )
%AK_GABAinASD_EyeTrackingFiguresGUI is a guided user interface designed to 
%read in specifically formatted data from a .xls file and turn it into a 
%beeswarm plot. Specifically, the function creates a subplot for each
%statistic selected via the GUI so that each condition of each scan type/ 
%experiment selected receives its own distribution spaced out along the x
%axis. Within each condition, there is the option to further separate data
%into groups (ASD vs NT for GABAinASD data; drug vs placebo for Lorazepam)
%so that each group * condition * experiment combination has a distribution 
%plotted on a common axis. There are also options to filter data to appear
%in different colors within each distribution based on subject ID or based
%on a criterion and some other statistic (i.e., the filter no tracking time
%> .8 would highlight points where no tracking time exceeds .8 in a
%different color than points which do not me that criterion). The GUI
%relies heavily on text labels that respond to mouse movements. Mouseover
%labels correspond to each individual data point, such that each point will
%be labeled as you hover a cursor over it. Also, while conventional
%x-labels describe subject groupings, sub-labels will pop up below the
%x-axis to describe what experiment and condition each distribution belongs
%to.
%
% Run the script Create_fMRI_EyeTracking_Data_Table_GABAinASD.m in order to
% generate correctly formatted .xls files for the GUI to read. The existing
% file(s) on the L drive should be named fMRI_EyeTracking_Data_Tables*.xlsx,
% fMRI_Behavioral_Data_Tables*.xlsx, or
% Psychophysics_EyeTracking_Data_Tables*.xlsx. Head motion data has its own
% simpler version of the GUI.
%   
%INPUT:
%   xlsfile [optional]: full file directory for .xls file where eye tracking data
%       tables are stored (string); defaults to UI selection from
%       directory: 'L:\MurrayLab\DataTablesForGUI'
%OUTPUT: (option to save output using 'Export' button)
%   exportData: a maxrix with 6 dimensions
%       (experiment,condition,statisitic,subject,session,set/run) or 7
%       dimensions
%       (experiment,condition,statisitic,subject,session,set/run,block),
%       the content of which will match whatever data set is currently
%       selected for visualization; blank indices are filled with nans
%   key: a structure containing descriptions of the data at each element of
%       exportData, listed in order of dimensions


%% check input 

if nargin < 1
    xlsfile_dir = 'L:\MurrayLab\DataTablesForGUI'; % set root dir (should be the same for any cpu where network drive is mapped to 'L:\')
    cd(xlsfile_dir);
    [xlsfile_name,xlsfile_dir] = uigetfile('*.xlsx','Select a file from which to load data');
    xlsfile = fullfile(xlsfile_dir,xlsfile_name); 
end

addpath(genpath('L:\MurrayLab\ASD\Data\Analysis_scripts\AMK_Code'));

%% initialize display and set up GUI referencing callback functions

% set prior values for variables spanning multiple functions
sheetsCurrent = {''}; % to pass between callback functions
condsCurrent = {''};
statsCurrent = {''};
subjectsCurrent = {''};
filterIdx = {[]};
filtstatSelected = '';
operatorSelected = '';
criterionSelected = [];
highlight_subjects = {''};
plotStruct = struct;

dataTables = {[]}; % to be defined in LoadData
sheetType = [];
condsAll = {''};
condsAllIdx = {[]};
condsAllWhichSheetIdx = [];
statsAll = {''};
statsAllIdx = {[]};
subjectsAll = {''};

plotByGroups = []; % boolean: whether or not to plot data by groups
exportData = []; % optional output (needs to have scope of full function)

% for reference
fMRIsessions = {'fMRI1','fMRI2'};
fMRIsets = {'set1','set2','set3','set4','set5','set6','set7','set8','set9'};
LZsessions = {'drug','placebo'};
PsychophysicsRuns = {'run1','run2','run3','run4'};

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
        uicontrol('tag','Experiments','string','Experiments','callback',@Spreadsheets_CB,'pos',[5 125 70 20]);
        uicontrol('tag','Statistics','string','Statistics','callback',@Statistics_CB,'pos',[5 102 60 20]);
        uicontrol('tag','Conditions','string','Conditions','callback',@Conditions_CB,'pos',[5 79 60 20]);
        uicontrol('tag','Subjects','string','Subjects','callback',@Subjects_CB,'pos',[5 56 60 20]);
        uicontrol('tag','Filter','string','Filter Data','callback',@Filter_CB,'pos',[5 33 60 20]);
        uicontrol('tag','ExportButton','string','Export','callback',@ExportData_CB,'pos',[5 10 60 20]);
        uicontrol('Style','togglebutton','tag','PlotByGroups','string','Plot by group','callback',@PlotByGroup_CB,'pos',[5 148 70 20])
    end
    
    function PlotByGroup_CB(~,~)
        % button compressed -> true; not compressed -> false
        buttonState = logical(get(findobj('tag','PlotByGroups'),'Value'));
        plotByGroups = buttonState;
        RePlot();
    end

%FYI: in order to set default data grouping, use the following code:
    % set(findobj('tag', 'PlotByGroups'), 'Value', [1 or 0])
    % feval(get(findobj('tag', 'PlotByGroups'), 'Callback'))

    function ExportData_CB(~,~)
        % get user input for where to save file
        [fileName,pathName] = uiputfile('*.mat', 'Save Current Selections as');
        % save formatted data
        if fileName ~= 0 % did they enter a name and directory?
            generalKey = {'experiments','conditions','statistics','subjects','sessions','set/runs','blocks'};
            experiment(1:length(sheetsCurrent)) = struct;
            for iExp = 1:length(sheetsCurrent);
                switch sheetType{iExp}
                    case {'block ftap', 'non-block ftap'}
                        experiment(iExp).experimentName = sheetsCurrent{iExp};
                        experiment(iExp).conditions = condsCurrent;
                        experiment(iExp).statistics = statsCurrent;
                        experiment(iExp).subjects = subjectsCurrent;
                        experiment(iExp).sessions = 'ftap';
                        experiment(iExp).sets = fMRIsets{1:3};
                        experiment(iExp).blocks = 'Blocks are in order from 1 to the maximum number of blocks run in a set. If this is a singleton dimension, there are no block-level statistics selected.';

                    case {'block fMRI', 'non-block fMRI'}
                        experiment(iExp).experimentName = sheetsCurrent{iExp};
                        experiment(iExp).conditions = condsCurrent;
                        experiment(iExp).statistics = statsCurrent;
                        experiment(iExp).subjects = subjectsCurrent;
                        experiment(iExp).sessions = fMRIsessions;
                        experiment(iExp).sets = fMRIsets;
                        experiment(iExp).blocks = 'Blocks are in order from 1 to the maximum number of blocks run in a set. If this is a singleton dimension, there are no block-level statistics selected.';
                    case 'GABA psychophysics'
                        experiment(iExp).experimentName = sheetsCurrent{iExp};
                        experiment(iExp).conditions = condsCurrent;
                        experiment(iExp).statistics = statsCurrent;
                        experiment(iExp).subjects = subjectsCurrent;
                        experiment(iExp).sessions = 'Sessions should be a singleton dimension for GABA psychophysics. All statistics for this experiment will be at an index of 1 in this dimension.';
                        experiment(iExp).sets = PsychophysicsRuns;
                        experiment(iExp).blocks = 'Blocks should be a singleton dimension for GABA psychophysics. All statistics for this experiment will be at an index of 1 in this dimension.';
                    case 'Lorazepam psychophysics'
                        experiment(iExp).experimentName = sheetsCurrent{iExp};
                        experiment(iExp).conditions = condsCurrent;
                        experiment(iExp).statistics = statsCurrent;
                        experiment(iExp).subjects = subjectsCurrent;
                        experiment(iExp).sessions = LZsessions;
                        experiment(iExp).sets = PsychophysicsRuns;
                        experiment(iExp).blocks = 'Blocks should be a singleton dimension for GABA psychophysics. All statistics for this experiment will be at an index of 1 in this dimension.';
                end
            end
            key = struct('general',generalKey,'experiment',experiment);
%             % old key
%             key.d1 = subjectsCurrent;
%             key.d2 = fMRIsessions;
%             key.d3 = {'set1','set2','set3','set4','set5','set6','set7','set8','set9'};
%             key.d4 = 'blocks';
%             key.d5 = condsCurrent;
%             key.d6 = statsCurrent;
            save(fullfile(pathName,fileName),'exportData','key')
        end
    end
    
    function LoadData() 
        % load eye tracking data tables
        for iSheet = 1:length(sheetsCurrent)
            % read in data
            [~,~,dataTables{iSheet}] = xlsread(xlsfile,sheetsCurrent{iSheet}); % read sheet
            % what kind of sheet is this? 
            % for the purpose of creating mouseover labels
            if strcmp(sheetsCurrent{iSheet},'block ftap') % is this a sheet for ftap block statistics?
                % label sheet type
                sheetType{iSheet} = 'block ftap';
            elseif strcmp(sheetsCurrent{iSheet},'ftap') % is this a sheet for non-block ftap statistics?
                % label sheet type
                sheetType{iSheet} = 'non-block ftap';
            elseif ~isempty(regexp(sheetsCurrent{iSheet},regexptranslate('wildcard','block*'),'ONCE')) % is this a sheet for fMRI block statistics?
                % label sheet type
                sheetType{iSheet} = 'block fMRI';
            elseif ~isempty(regexp(sheetsCurrent{iSheet},regexptranslate('wildcard','*GABA*'),'ONCE')) % is this a sheet for GABA psychophysics statistics?
                % label sheet type
                sheetType{iSheet} = 'GABA psychophysics';
            elseif ~isempty(regexp(sheetsCurrent{iSheet},regexptranslate('wildcard','*LZ*'),'ONCE')) % is this a sheet for Lorazepam psychophysics statistics?
                % label sheet type
                sheetType{iSheet} = 'Lorazepam psychophysics';
            else
                % label sheet type
                sheetType{iSheet} = 'non-block fMRI';
            end
            % fill in empty cells with nans
            dataTables{iSheet}(cellfun(@isempty,dataTables{iSheet})) = {nan};
            
            % find subjects list
            tableSubjects{iSheet} = unique(dataTables{iSheet}(2:end,1));
            % find condition and stat names for current sheet
            clear condStrIndex tableCondListFull tableStatListFull tableCondListUnique condIdx statIdx tableCondList
            findCondStatStr = '*:*';
            condStatStrIndex = find(~cellfun(@isempty,regexp(dataTables{iSheet}(1,:),regexptranslate('wildcard',findCondStatStr)))==1);
            for iCS = 1:length(condStatStrIndex)
                clear tempIdx
                tempIdx = strfind(dataTables{iSheet}{1,condStatStrIndex(iCS)},findCondStatStr(2)); % idx for position of : in each string where : occurs
                tableCondListFull{iCS} = dataTables{iSheet}{1,condStatStrIndex(iCS)}(1:tempIdx-1); % list of all condition names in order of columns
                tableStatListFull{iCS} = dataTables{iSheet}{1,condStatStrIndex(iCS)}(tempIdx+2:end); % list of all stat names in order of columns
            end
            tableCondListUnique = unique(tableCondListFull,'stable'); % unique cond names
            tableStatList{iSheet} = unique(tableStatListFull,'stable'); % unique stat names
            % create column index for each condition on current sheet
            for iC = 1:length(tableCondListUnique);
                condIdx{iC} = find(strcmp(tableCondListFull,tableCondListUnique{iC})) + condStatStrIndex(1) - 1; % create an index for columns for each condition
            end
            % create column index for each stat on current sheet
            for iS = 1:length(tableStatList{iSheet});
                statIdx{iS} = find(~cellfun(@isempty,regexp(dataTables{iSheet}(1,:),regexptranslate('wildcard',['*: ' tableStatList{iSheet}{iS}])))); % create an index for columns for each stat
            end
            % add conds and indices for current table to condsAll and condsAllIdx lists
            tableCondList = cellfun(@(c) [sheetsCurrent{iSheet} ':' c],tableCondListUnique,'UniformOutput',false); % label each condition with its spreadsheet name
            if all(cellfun(@isempty,condsAll))
                condsAll = tableCondList';
            else
                condsAll = [condsAll; tableCondList'];
            end
            if all(cellfun(@isempty,condsAllIdx))
                condsAllIdx = condIdx';
            else
                condsAllIdx = [condsAllIdx; condIdx'];
            end
            if isempty(condsAllWhichSheetIdx)
                condsAllWhichSheetIdx = repmat(iSheet,[length(tableCondList) 1]);
            else
                condsAllWhichSheetIdx = [condsAllWhichSheetIdx; repmat(iSheet,[length(tableCondList) 1])];
            end
            % add stat indices for current table to statsAllIdx list
            statsAllIdx{iSheet} = statIdx;
        end
        % store subject codes and stat names subjectsAll and statsAll list
        if iSheet>1
            if isequal(tableSubjects{1:iSheet}) % use first list if all are identical
                subjectsAll = tableSubjects{1};
            else % if not, use the shorter list of subjects
                warning('The selected sheets do not contain all the same subjects.')
                uselistIdx = cellfun(@length,tableSubjects)==min(cellfun(@length,tableSubjects));
                subjectsAll = tableSubjects{uselistIdx};
                % chech that there are no subjects in list that are not in all sheets
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
        condsCurrent = condsAll;
    end

    function Spreadsheets_CB(~,~)
        % select spreadsheets
        [~,sheetsAll] = xlsfinfo(xlsfile);
        sheetsAll(strcmp('Sheet1',sheetsAll)) = []; % remove Sheet1 (excel default and an empty sheet) as an option
        if all(cellfun(@isempty,sheetsCurrent))
            sheetsSelected = sheetsAll(listdlg('liststring',sheetsAll,'PromptString','Select Experiment(s)'));
            sheetsCurrent = sheetsSelected;
            % reset other variables to prior values
            condsCurrent = {''}; % to pass between callback functions
            statsCurrent = {''};
            subjectsCurrent = {''};
            filterIdx = {[]};
            filtstatSelected = '';
            operatorSelected = '';
            criterionSelected = [];
            highlight_subjects = {''};
            dataTables = {[]}; % to be defined in LoadData
            sheetType = [];
            condsAll = {''};
            condsAllIdx = {[]};
            condsAllWhichSheetIdx = [];
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
                condsCurrent = {''}; % to pass between callback functions
                statsCurrent = {''};
                subjectsCurrent = {''};
                filterIdx = {''};
                filtstatSelected = '';
                operatorSelected = '';
                criterionSelected = [];
                highlight_subjects = {''};
                dataTables = {[]}; % to be defined in LoadData
                sheetType = [];
                condsAll = {''};
                condsAllIdx = {[]};
                condsAllWhichSheetIdx = [];
                statsAll = {''};
                statsAllIdx = {[]};
                subjectsAll = {''};
                if ~all(cellfun(@isempty,sheetsCurrent))
                    LoadData();
                end
            end
        end
    end

    function Conditions_CB(~,~)
        % must run spreadsheets if condsAll is empty
        if all(cellfun(@isempty,condsAll))
            Spreadsheets_CB();
        end
        % select conditions
        if all(cellfun(@isempty,condsCurrent))
            condsSelected = condsAll(listdlg('liststring',condsAll,'PromptString','Select Condition(s)'));
            condsCurrent = condsSelected; 
            if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                RePlot(); 
            end
        else
            condsSelected = condsAll(listdlg('liststring',condsAll,'initialvalue',find(ismember(condsAll,condsCurrent)),'PromptString','Select Condition(s)'));
            % replot?
            if numel(condsCurrent) ~= numel(condsSelected) || ~all(ismember(condsCurrent,condsSelected)) || ~all(ismember(condsSelected,condsCurrent))
               condsCurrent = condsSelected; 
               if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                    RePlot(); 
                end 
            end
        end
    end

    function Statistics_CB(~,~)
         % must run spreadsheets if statsAll is empty
        if all(cellfun(@isempty,statsAll))
            Spreadsheets_CB();
        end
        % select statistics
        if all(cellfun(@isempty,statsCurrent))
            statsSelected = statsAll(listdlg('liststring',statsAll,'PromptString','Select Statistic(s)'));
            statsCurrent = statsSelected;
            if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                RePlot(); 
            end 
        else
            statsSelected = statsAll(listdlg('liststring',statsAll,'initialvalue',find(ismember(statsAll,statsCurrent)),'PromptString','Select Statistic(s)'));
            % replot?
            if numel(statsCurrent) ~= numel(statsSelected) || ~all(ismember(statsCurrent,statsSelected)) || ~all(ismember(statsSelected,statsCurrent))
               statsCurrent = statsSelected; 
               if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
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
            if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                RePlot(); 
            end 
        else
            subjectsSelected = subjectsAll(listdlg('liststring',subjectsAll,'initialvalue',find(ismember(subjectsAll,subjectsCurrent)),'PromptString','Select Subject(s)'));
            % replot?
            if numel(subjectsCurrent) ~= numel(subjectsSelected) || ~all(ismember(subjectsCurrent,subjectsSelected)) || ~all(ismember(subjectsSelected,subjectsCurrent))
               subjectsCurrent = subjectsSelected;
               highlight_subjects = {[]}; % reset highlighted subjects
               if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
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
            if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
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
            if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
                RePlot(); 
            end 
        else
            subjectsSelected = subjectsCurrent(listdlg('liststring',subjectsCurrent,'initialvalue',find(ismember(subjectsCurrent,highlight_subjects)),'PromptString','Select Subject(s) to Highlight'));
            % replot?
            if numel(highlight_subjects) ~= numel(subjectsSelected) || ~all(ismember(highlight_subjects,subjectsSelected)) || ~all(ismember(subjectsSelected,highlight_subjects))
               highlight_subjects = subjectsSelected; 
               if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
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
                statColumnIdx = ~cellfun(@isempty,regexp(dataTables{iSheet}(1,:),regexptranslate('wildcard',['*: ' filtstatSelected])));
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
                if ~all(cellfun(@isempty,statsCurrent)) && ~all(cellfun(@isempty,condsCurrent)) && ~all(cellfun(@isempty,subjectsCurrent))
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

% change the way that groups are indexed to be more flexible:
    % should accomodate any number of groups greater than 1
    % should depend on value of plotByGroups and sheetType variables
    

    function RePlot()
        % separate into vectors of data by group, stats of interest, generate labels, and group information for plotting
        for iSheet = 1:length(sheetsCurrent);
            % indices for group membership: groupIdx{sheet/experiment, group}
            if plotByGroups
                % parse groups depending on sheet type
                switch sheetType{iSheet}
                    case 'Lorazepam psychophysics'
                        groupIdx{iSheet,1} = find(strcmp(dataTables{iSheet}(:,2),'drug')); % group1 == drug
                        groupIdx{iSheet,2} = find(strcmp(dataTables{iSheet}(:,2),'placebo')); % group2 == placebo
                        groupLabels{iSheet} = {'drug','placebo'}; % to use in x-axis labels
                    otherwise % GABA in ASD 
                        groupIdx{iSheet,1} = find(~cellfun(@isempty,regexp(dataTables{iSheet}(:,1),regexptranslate('wildcard','*G1*')))==1); % group1 == ASD
                        groupIdx{iSheet,2} = find(~cellfun(@isempty,regexp(dataTables{iSheet}(:,1),regexptranslate('wildcard','*G3*')))==1); % group2 == NT
                        groupLabels{iSheet} = {'ASD','NT'}; % to use in x-axis labels
                end
            else
                % don't index by groups, but maintain similar data indexing
                groupIdx{iSheet,1} = 2:length(dataTables{iSheet}(:,1)); % group1 == all subjects; exclude first row of column labels from index
                groupLabels{iSheet} = {'All'}; % to use in x-axis labels
            end
            % use only subjects currently selected
            if ~all(cellfun(@isempty,subjectsCurrent))
                for iGroup = 1:length(groupIdx(iSheet,:))
                    groupIdx{iSheet,iGroup} = groupIdx{iSheet,iGroup}(ismember(dataTables{iSheet}(groupIdx{iSheet,iGroup},1),subjectsCurrent));
                end
            end
            % index for subjects of interest (highlight_subjects)
            if ~all(cellfun(@isempty,highlight_subjects))
                for iGroup = 1:length(groupIdx(iSheet,:))
                    highlightIdx{iSheet,iGroup} = zeros(length(groupIdx{iSheet,iGroup}),1);
                    for iSubj = 1:length(highlight_subjects);
                        highlightIdx{iSheet,iGroup} = highlightIdx{iSheet,iGroup} + strcmp(dataTables{iSheet}(groupIdx{iSheet,iGroup},1),highlight_subjects{iSubj});
                    end
                    highlightIdx{iSheet,iGroup} = logical(highlightIdx{iSheet,iGroup});
                end
            end
        end
        
        % preallocate counter
        x_labelIdx = 1;
        
        % assignment of data arrays
        for iCond = 1:length(condsCurrent)
            % index for this cond and the associated sheet
            clear thisCondIdx thisCondSheetIdx
            thisCondIdx = strcmp(condsAll,condsCurrent{iCond}); % idx in condsAll associated with iCond
            thisCondSheetIdx = condsAllWhichSheetIdx(strcmp(condsAll,condsCurrent{iCond})); % sheet idx associated with iCond
            for iStat = 1:length(statsCurrent)
                % index for this stat 
                clear otherStats thisStatIdx
                otherStats = setdiff(statsAll,statsCurrent{iStat});
                thisStatIdx = strcmp(statsAll,statsCurrent{iStat}) & ~cellfun(@(x) any(strcmp(otherStats,x)),statsAll); % idx in statsAll associated with iStat
                % store column values by stat and cond
                clear statCondIdx
                statCondIdx = intersect(statsAllIdx{thisCondSheetIdx}{thisStatIdx},condsAllIdx{thisCondIdx}); % should be only one column where stat and cond intersect
                % assign data to plot to stat structure
                stat(iStat).name = statsCurrent{iStat};
                for iGroup = 1:length(groupIdx(thisCondSheetIdx,:))
                    stat(iStat).group(iGroup).values{iCond} = cell2mat(dataTables{thisCondSheetIdx}(groupIdx{thisCondSheetIdx,iGroup},statCondIdx));
                end
            
                % create mouseover labels for subject,fMRI session, and
                % set# and prepare exportData matrix
                switch sheetType{thisCondSheetIdx} % labels different for block and summary stats
                    case 'non-block ftap' % for summary stats:
                        % mouseover labels and temporary arrays
                        tempLabels = cell(0); % preallocate
                        tempData = [];
                        for iGroup = 1:length(groupIdx(thisCondSheetIdx,:))
                            % mouseover labels
                            for iLabel = 1:length(groupIdx{thisCondSheetIdx,iGroup})
                                stat(iStat).group(iGroup).labels{iCond}{iLabel} = [dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),1} ': ' dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),2} ': ' dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),3}];
                            end
                            % temporary arrays for data export
                            tempLabels = [tempLabels stat(iStat).group(iGroup).labels{iCond}];
                            tempData = [tempData; stat(iStat).group(iGroup).values{iCond}];
                        end
                        
                        %%% format data for export (01/26/17):
                        for iSubj = 1:length(subjectsCurrent)
                            for iSess = 1
                                for iSet = 1:3
                                    % indices for this subject, session, and set
                                    if isfield(stat(iStat).group,'labels')
                                        thisSubjIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
                                        thisSessIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard','*ftap*')));
                                        thisSetIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*set' num2str(iSet) '*'])));
                                    else
                                        warning('No subjects of either group selected');
                                        thisSubjIdx = false;
                                        thisSessIdx = false;
                                        thisSetIdx = false;
                                    end
                                    % store data
                                    if isempty(find(thisSubjIdx & thisSessIdx & thisSetIdx,1))
                                        exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iSet,1) = nan;
                                    else
                                        exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iSet,1) = tempData(thisSubjIdx & thisSessIdx & thisSetIdx);
                                    end
                                end
                            end
                        end
                    case 'block ftap' % for block stats:
                        % mouseover labels and temporary arrays
                        tempLabels = cell(0); % preallocate
                        tempData = [];
                        maxBlocks = 0;
                        for iGroup = 1:length(groupIdx(thisCondSheetIdx,:))
                            for iLabel = 1:length(groupIdx{thisCondSheetIdx,iGroup})
                                % mouseover labels for plots
                                stat(iStat).group(iGroup).labels{iCond}{iLabel} = [dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),1} ': ' dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),2} ': ' dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),3} ': block' num2str(dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),5})];
                                % for export
                                if maxBlocks < dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),5}
                                    maxBlocks = dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),5};
                                end
                                % temporary arrays for data export
                                tempLabels = [tempLabels stat(iStat).group(iGroup).labels{iCond}];
                                tempData = [tempData; stat(iStat).group(iGroup).values{iCond}];
                            end
                        end
                        
                        %%% format data for export (01/26/17):
                        for iSubj = 1:length(subjectsCurrent)
                            for iSess = 1
                                for iSet = 1:3
                                    for iBlock = 1:maxBlocks
                                        % indices for this subject, session, and set
                                        if isfield(stat(iStat).group,'labels')
                                            thisSubjIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
                                            thisSessIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard','*ftap*')));
                                            thisSetIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*set' num2str(iSet) '*'])));
                                            thisBlockIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*block' num2str(iBlock) '*'])));
                                        else
                                            warning('No subjects of either group selected');
                                            thisSubjIdx = false;
                                            thisSessIdx = false;
                                            thisSetIdx = false;
                                            thisBlockIdx = false;
                                        end
                                        % store data
                                        if isempty(find(thisSubjIdx & thisSessIdx & thisSetIdx & thisBlockIdx,1))
                                            exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iSet,iBlock) = nan;
                                        else
                                            exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iSet,iBlock) = tempData(thisSubjIdx & thisSessIdx & thisSetIdx & thisBlockIdx);
                                        end
                                    end
                                end
                            end
                        end
                    case 'non-block fMRI' % for summary stats:
                        % mouseover labels and temporary arrays
                        tempLabels = cell(0); % preallocate
                        tempData = [];
                        for iGroup = 1:length(groupIdx(thisCondSheetIdx,:))
                            for iLabel = 1:length(groupIdx{thisCondSheetIdx,iGroup})
                                % mouseover labels
                                stat(iStat).group(iGroup).labels{iCond}{iLabel} = [dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),1} ': ' dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),2} ': ' dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),3}];
                                % temporary arrays for data export
                                tempLabels = [tempLabels stat(iStat).group(iGroup).labels{iCond}];
                                tempData = [tempData; stat(iStat).group(iGroup).values{iCond}];
                            end
                        end
                        
                        %%% format data for export (01/26/17):
                        for iSubj = 1:length(subjectsCurrent)
                            for iSess = 1:2
                                for iSet = 1:9
                                    % indices for this subject, session, and set
                                    if isfield(stat(iStat).group,'labels')
                                        thisSubjIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
                                        thisSessIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*fMRI' num2str(iSess) '*'])));
                                        thisSetIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*set' num2str(iSet) '*'])));
                                    else
                                        warning('No subjects of either group selected');
                                        thisSubjIdx = false;
                                        thisSessIdx = false;
                                        thisSetIdx = false;
                                    end
                                    % store data
                                    if isempty(find(thisSubjIdx & thisSessIdx & thisSetIdx,1))
                                        exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iSet,1) = nan;
                                    else
                                        exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iSet,1) = tempData(thisSubjIdx & thisSessIdx & thisSetIdx);
                                    end
                                end
                            end
                        end
                    case 'block fMRI' % for block stats:
                        % mouseover labels and temporary arrays
                        tempLabels = cell(0); % preallocate
                        tempData = [];
                        maxBlocks = 0; 
                        for iGroup = 1:length(groupIdx(thisCondSheetIdx,:))
                            for iLabel = 1:length(groupIdx{thisCondSheetIdx,iGroup})
                                % mouseover labels for plots
                                stat(iStat).group(iGroup).labels{iCond}{iLabel} = [dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),1} ': ' dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),2} ': ' dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),3} ': block' num2str(dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),5})];
                                % for export
                                if maxBlocks < dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),5}
                                    maxBlocks = dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),5};
                                end
                                % temporary arrays for data export
                                tempLabels = [tempLabels stat(iStat).group(iGroup).labels{iCond}];
                                tempData = [tempData; stat(iStat).group(iGroup).values{iCond}];
                            end
                        end
                        
                        %%% format data for export (01/26/17):
                        for iSubj = 1:length(subjectsCurrent)
                            for iSess = 1:2
                                for iSet = 1:9
                                    for iBlock = 1:maxBlocks
                                        % indices for this subject, session, and set
                                        if isfield(stat(iStat).group,'labels')
                                            thisSubjIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
                                            thisSessIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*fMRI' num2str(iSess) '*'])));
                                            thisSetIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*set' num2str(iSet) '*'])));
                                            thisBlockIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*block' num2str(iBlock) '*'])));
                                        else
                                            warning('No subjects of either group selected');
                                            thisSubjIdx = false;
                                            thisSessIdx = false;
                                            thisSetIdx = false;
                                            thisBlockIdx = false;
                                        end
                                        % store data
                                        if isempty(find(thisSubjIdx & thisSessIdx & thisSetIdx & thisBlockIdx,1))
                                            exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iSet,iBlock) = nan;
                                        else
                                            exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iSet,iBlock) = tempData(thisSubjIdx & thisSessIdx & thisSetIdx & thisBlockIdx);
                                        end
                                    end
                                end
                            end
                        end
                    case 'GABA psychophysics'
                        % mouseover labels and temporary arrays
                        tempLabels = cell(0); % preallocate
                        tempData = [];
                        for iGroup = 1:length(groupIdx(thisCondSheetIdx,:))
                            % mouseover labels
                            for iLabel = 1:length(groupIdx{thisCondSheetIdx,iGroup})
                                stat(iStat).group(iGroup).labels{iCond}{iLabel} = [dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),1} ': run' num2str(dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),2})];
                            end
                            % temporary arrays for data export
                            tempLabels = [tempLabels stat(iStat).group(iGroup).labels{iCond}];
                            tempData = [tempData; stat(iStat).group(iGroup).values{iCond}];
                        end
                        
                        %%% format data for export (01/26/17):
                        for iSubj = 1:length(subjectsCurrent)
                            for iRun = 1:4
                                % indices for this subject, session, and set
                                if isfield(stat(iStat).group,'labels')
                                    thisSubjIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
                                    thisRunIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*run' num2str(iRun) '*'])));
                                else
                                    warning('No subjects of either group selected');
                                    thisSubjIdx = false;
                                    thisRunIdx = false;
                                end
                                % store data
                                if isempty(find(thisSubjIdx & thisRunIdx,1))
                                    exportData(thisCondSheetIdx,iCond,iStat,iSubj,1,iRun,1) = nan;
                                else
                                    exportData(thisCondSheetIdx,iCond,iStat,iSubj,1,iRun,1) = tempData(thisSubjIdx & thisRunIdx);
                                end
                            end
                        end
                    case 'Lorazepam psychophysics'
                        % mouseover labels and temporary arrays
                        tempLabels = cell(0); % preallocate
                        tempData = [];
                        for iGroup = 1:length(groupIdx(thisCondSheetIdx,:))
                            % mouseover labels
                            for iLabel = 1:length(groupIdx{thisCondSheetIdx,iGroup})
                                stat(iStat).group(iGroup).labels{iCond}{iLabel} = [dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),1} ': ' dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),2} ': run' num2str(dataTables{thisCondSheetIdx}{groupIdx{thisCondSheetIdx,iGroup}(iLabel),3})];
                            end
                            % temporary arrays for data export
                            tempLabels = [tempLabels stat(iStat).group(iGroup).labels{iCond}];
                            tempData = [tempData; stat(iStat).group(iGroup).values{iCond}];
                        end
                        
                        %%% format data for export (01/26/17):
                        for iSubj = 1:length(subjectsCurrent)
                            for iSess = 1:2 % drug (1) vs placebo (2)
                                for iRun = 1:4
                                    % indices for this subject, session, and set
                                    if isfield(stat(iStat).group,'labels')
                                        thisSubjIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*' subjectsCurrent{iSubj} '*'])));
                                        thisSessIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*' LZsessions{iSess} '*'])));
                                        thisRunIdx = ~cellfun(@isempty,regexp(tempLabels,regexptranslate('wildcard',['*run' num2str(iRun) '*'])));
                                    else
                                        warning('No subjects of either group selected');
                                        thisSubjIdx = false;
                                        thisSessIdx = false;
                                        thisRunIdx = false;
                                    end
                                    % store data
                                    if isempty(find(thisSubjIdx & thisSessIdx & thisRunIdx,1))
                                        exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iRun,1) = nan;
                                    else
                                        exportData(thisCondSheetIdx,iCond,iStat,iSubj,iSess,iRun,1) = tempData(thisSubjIdx & thisSessIdx & thisRunIdx);
                                    end
                                end
                            end
                        end
                end
            end
            
            % use filters only for selected conds 
            if ~all(cellfun(@isempty,filterIdx))
                % generate index for the condition number wihtin the current sheet for accessing filterIdx (easier to do this here rather than later)
                clear thisSheetCondIdx
                thisSheetCondIdx = strcmp(condsAll(~cellfun(@isempty,regexp(condsAll,regexptranslate('wildcard',[sheetsCurrent{thisCondSheetIdx} ':'])))),condsCurrent{iCond});
                % set up useFiltIdx as a cell array with one cell per iCond
                % this will allow for easier indexing in the creation of highlightIndices
                useFiltIdx{iCond} = filterIdx{thisCondSheetIdx}(:,thisSheetCondIdx);
            end

            % x-axis labels:
            % do this differently depending on how selected subjects are grouped
            for iGroup = 1:length(groupIdx(thisCondSheetIdx,:))
                if any(cellfun(@any,groupIdx(thisCondSheetIdx,iGroup)))
                    x_labels{x_labelIdx} = [condsCurrent{iCond} ' ' groupLabels{thisCondSheetIdx}{iGroup}];
                    x_labelIdx = x_labelIdx + 1;
                end
            end
            x_labels{x_labelIdx} = ''; % blank label for spacing between coonditions
            x_labelIdx = x_labelIdx + 1;
        end
        x_labels(end) = []; % remove last blank x-axis label
        
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
                   clear labelCondIdx labelSheetIdx labelGroupIdx
                   % do this differently for block and non-block x_labels
                   if ~isempty(strfind(x_labels{iX},'block')) % if this x_label is for stats by block
                       labelCondIdx = find(AK_whichPattern(x_labels{iX},condsCurrent) & ~cellfun(@isempty,strfind(condsCurrent,'block'))'); % which cond corresponds to x_labels{iX}
                       labelSheetIdx = find(AK_whichPattern(x_labels{iX},sheetsCurrent) & ~cellfun(@isempty,strfind(sheetsCurrent,'block'))); % which sheet corresponds to x_labels{iX}
                   else % if this x_label is for stats by set, rather than by block
                       labelCondIdx = find(AK_whichPattern(x_labels{iX},condsCurrent) & cellfun(@isempty,strfind(condsCurrent,'block'))'); % which cond corresponds to x_labels{iX}
                       labelSheetIdx = find(AK_whichPattern(x_labels{iX},sheetsCurrent) & cellfun(@isempty,strfind(sheetsCurrent,'block'))); % which sheet corresponds to x_labels{iX}
                   end
                   % what group does this data belong to?
                   labelGroupIdx = find(AK_whichPattern(x_labels{iX},groupLabels{labelSheetIdx}));
                   % does this label match any group and is there any data for this intersection of stat, group, and condition?
                   if ~isempty(labelGroupIdx) && ~cellfun(@isempty,stat(iStat).group(labelGroupIdx).values(labelCondIdx)) 
                       % match group x_labels for cond with proper data arrays and empty datapoint labels
                       stat(iStat).plotData(iX) = stat(iStat).group(labelGroupIdx).values(labelCondIdx);
                       stat(iStat).plotLabels(iX) = stat(iStat).group(labelGroupIdx).labels(labelCondIdx);
                       stat(iStat).x_labels(iX) = groupLabels{labelSheetIdx}(labelGroupIdx);
                       % append highlightIndices
                       if ~all(cellfun(@isempty,highlight_subjects)) && all(cellfun(@isempty,filterIdx)) && iStat==1
                           highlightIndices = vertcat(highlightIndices,highlightIdx{labelSheetIdx,labelGroupIdx});
                           distIdx = vertcat(distIdx,repmat(iX,size(highlightIdx{labelSheetIdx,labelGroupIdx})));
                       end
                       if all(cellfun(@isempty,highlight_subjects)) && ~all(cellfun(@isempty,filterIdx)) && iStat==1
                           highlightIndices = vertcat(highlightIndices,2.*useFiltIdx{labelCondIdx}(groupIdx{labelSheetIdx,labelGroupIdx}));
                           distIdx = vertcat(distIdx,repmat(iX,size(2.*useFiltIdx{labelCondIdx}(groupIdx{labelSheetIdx,labelGroupIdx}))));
                       end
                       if ~all(cellfun(@isempty,highlight_subjects)) && ~all(cellfun(@isempty,filterIdx)) && iStat==1
                           highlightIndices = vertcat(highlightIndices, highlightIdx{labelSheetIdx,labelGroupIdx} + 2.*useFiltIdx{labelCondIdx}(groupIdx{labelSheetIdx,labelGroupIdx}) );
                           distIdx = vertcat(distIdx,repmat(iX,size(highlightIdx{labelSheetIdx,labelGroupIdx} + 2.*useFiltIdx{labelCondIdx}(groupIdx{labelSheetIdx,labelGroupIdx}))));
                       end
                       % append x axis label
                       if iStat==1
                           xAxisLabel = [xAxisLabel {AK_eraseSubstring(x_labels{iX},[' ' groupLabels{labelSheetIdx}{labelGroupIdx}])}]; % add experiment name to x axis label
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
            % hack to make sure legend is correct
%             subplot(length(statsCurrent),1,1);
            % need to plot elements in order and in a different quadrant of
            % the axes than the rest of the data and then set the axes
            % limits to the first quadrant
            for iStat = 1:length(statsCurrent)
                % plot
                subplot(length(statsCurrent),1,iStat);
%                 stat(iStat).handles = AK_plotSpread_mouseover(stat(iStat).plotData,stat(iStat).plotLabels,'categoryIdx',highlightIndices,'xNames',stat(iStat).x_labels,'showMM',5);
                stat(iStat).handles = AK_plotSpread_forGUI(stat(iStat).plotData,'categoryIdx',highlightIndices,'xNames',stat(iStat).x_labels,'showMM',5);                
                title(stat(iStat).name);
                set(stat(iStat).handles{2}(2),'LineWidth',1.5)
                % figure out y position of x labels relative to axes
                yL = get(gca,'YLim');
                xAxisLabelY = yL(1) - (0.1 + .07*(length(statsCurrent)-1))*(yL(2)-yL(1)); % lower limit of y axis minus 10% of y axis range, minus an extra 7% per subplot
                % experiment labels as text
                for iE = 1:length(xAxisLabel)
                    t = text(xAxisLabelX(iE),xAxisLabelY,xAxisLabel{iE},'tag','xlabel experiment');
                    set(t,'Visible','off');
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
                set(stat(iStat).handles{2}(2),'LineWidth',1.5)
                % figure out y position of x labels relative to axes
                yL = get(gca,'YLim');
                xAxisLabelY = yL(1) - (0.1 + .07*(length(statsCurrent)-1))*(yL(2)-yL(1)); % lower limit of y axis minus 10% of y axis range, minus an extra 7% per subplot
                % experiment labels as text
                for iE = 1:length(xAxisLabel)
                    t = text(xAxisLabelX(iE),xAxisLabelY,xAxisLabel{iE},'tag','xlabel experiment');
                    set(t,'Visible','off');
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
            % pass input to single mousemove function that deals with both
            % mouseover labels and x-axis label visibility
            set(fH,'windowbuttonmotionfcn',{@mousemove,plotH,rawData,mouseoverLabels});
            % pass input to function which redraws x-labels when window is resized
            set(fH,'ResizeFcn',{@drawXLabels,xAxisLabel,xAxisLabelX})
        end   
    end
    
    function mousemove(src,ev,plotHandles,data,labels)
        %since this is a figure callback, the first input is the figure handle:
        f=handle(src);
        %like all callbacks, the second input, ev, isn't used.
        
        %determine which object is below the cursor:
        obj=hittest(f); %<-- the important line in this function
        
        % call mouseover labels function
        mouseoverLabels(f,obj,plotHandles,data,labels);
        % call xlabel visibility function
        xLabelVisibility(f,obj)
    end

    function mouseoverLabels(figHandle,currentObj,plotHandles,data,labels)
        if ~isequal(currentObj,figHandle) % if the current object is not the figure window
            if any(any(any(cellfun(@(x) isequal(currentObj,x),plotHandles)))) %if over any plot... 
                % change pointer when hovering over plot
                set(figHandle,'Pointer','crosshair');

                % determine which stat and distribution plot
                plotHindex = find(cellfun(@(x) isequal(currentObj,x),plotHandles));
                [~, plotHcol, ~] = ind2sub(size(plotHandles),plotHindex);

                %get cursor coordinates in its axes:
                a=get(plotHandles{plotHindex},'parent');
                point=get(a,'currentpoint');
                xclick=point(1,1,1);
                yclick=point(1,2,1);

                %determine which point we're over:
                idx=findclosestpoint2D(xclick,yclick,plotHandles{plotHindex});

                %make a "tool tip" that displays data point
                xl = xlim; yl = ylim;
                xrange = xl(2)-xl(1); yrange = yl(2)-yl(1);
                xoffset=xrange/1000;
                yoffset=yrange/1000;

                delete(findobj(figHandle,'tag','mouseover')); %delete last mouseover
                text(plotHcol+xoffset,data{plotHindex}(idx)+yoffset,labels{plotHindex}{idx},...
                    'backgroundcolor',[1 1 .8],'tag','mouseover','edgecolor',[0 0 0],...
                    'hittest','off');
            else
                % change pointer back and delete plot
                set(gcf,'Pointer','arrow');
                delete(findobj(figHandle,'tag','mouseover')); %delete last mouseover
            end
        else
            % change pointer back and delete plot
            set(gcf,'Pointer','arrow');
            delete(findobj(figHandle,'tag','mouseover')); %delete last mouseover
        end
    end

    function xLabelVisibility(figHandle,currentObj)
        % get current axes
        a = findobj(figHandle,'Type','Axes');
        if ~isequal(currentObj,figHandle) % if the current object is not the figure window
            % what is this object?
            objIDidx = cellfun(@(x) arrayfun(@(y) isequal(currentObj,y),x),... % for each set of axis children, return an array of logicals looking for a match for currentObj
                arrayfun(@allchild,a,'UniformOutput',false),... % returns a cell array where each cell contains an array of children corresponding to each axis in this figure 
                'UniformOutput',false);
            % is it a child of any axis on this figure?
            axisIdx = find(cellfun(@any, objIDidx));
            if ~isempty(axisIdx) %if over any axis?
                %get cursor coordinates in its axes:
                point=get(a(axisIdx),'currentpoint');
                xclick=point(1,1,1);
                yclick=point(1,2,1);

                % text handles for xlabel for experiments on the current axis
                t = findobj(a(axisIdx),'tag','xlabel experiment');
                
                % determine which x-label is closest
                idx = findclosestpoint2D(xclick,yclick,t);
                
                % make all x axis labels invisible
                for iA = 1:length(a)
                    allT = findobj(a(iA),'Type','Text');
                    % exclude current mouseover label from this list without
                    % eleminating mouseover label text object
                    mouseoverIdx = arrayfun(@(x) isequal(x,findobj(allT,'tag','mouseover')),allT);
                    copyT = allT(~mouseoverIdx);
                    for iE = 1:length(copyT)
                        set(copyT(iE),'Visible','off');
                    end
                end
                
                % make current x axis label visbile
                set(t(idx),'Visible','on');
            end
        end      
    end

    function index=findclosestpoint2D(xclick,yclick,datasource)
        %this function checks which point in the plotted line "datasource"
        %is closest to the point specified by xclick/yclick.

        if all(isprop(datasource,'xdata') & isprop(datasource,'ydata'))
            % for data points
            xdata = get(datasource,'xdata');
            ydata = get(datasource,'ydata');
        elseif all(isprop(datasource,'position'))
            % for non-data points
            positions = get(datasource,'position');
            xdata = cellfun(@(x) x(1), positions);
            ydata = cellfun(@(x) x(2), positions);
        end

        activegraph=get(datasource,'parent');
        
        if iscell(activegraph)
            % deal with the case where get parent returns a cell array of
            % references, all to the same parent axis
            activegraph = activegraph{1}; 
        end

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
            xAxisLabelY = yL(1) - (0.1 + .07*(length(a)-1))*(yL(2)-yL(1)); % lower limit of y axis minus 10% of y axis range, minus an extra 7% per subplot
            % text handles for xlabel for experiments on the current axis
            t = findobj(a(iA),'tag','xlabel experiment');
            % experiment labels as text
            for iE = 1:length(xAxisLabel)
                t(iE) = text(xAxisLabelX(iE),xAxisLabelY,xAxisLabel{iE},'tag','xlabel experiment');
                set(t(iE),'Visible','off');
            end
        end
    end

end

