function [ structure ] = AK_SupSum_makeHeadMotionTableFigures( sheetName, highlight_subjects, xlsfile_dir )
%AK_SupSum_makeEyetrackingTableFigures uses eye tracking data stored in 
%'L:\MurrayLab\ASD\Data\fMRI_HeadMotion_Data_Tables.xlsx'to create 
%beeswarm plots. Specifically, the function creates a figures for
%the type of scan denoted in 'sheetName', and it creates as many figures as
%the number of stats stored in the table. Within these figures, each group 
%(ASD and NT) has a distribution plotted on a common axis. There is also
%the option to highlight the data for specific subjects in green.
%   
%INPUT:
%   sheetName: the name of the excel sheet you wish to visualize (string);
%       must be an exact match
%   highlight_subjects [optional]: the name of the subjects you would like
%       to see highlighted in green on the beeswarm plots (string or cell array of
%       strings); defaults to none
%   xlsfile_dir[optional]: full file directory for xls file where eye tracking data
%       tables are stored (string); defaults to 'L:\MurrayLab\ASD\SuppressionSummation\fMRI_HeadMotion_Data_Tables.xlsx'
%OUTPUT:
%   structure = 
%       statName: same as 'plotStats'
%       handle: for each figure
%       plotData: the data plotted in the figures, separated by group and
%           condition
%       x_labels: x-axis labels (correspond to dimensions of 'plotData' and
%           'mouseover_labels)
%       mouseover_labels: labels associated with each plotted data point, 
%           separated by group and condition


%% check inputs 

if nargin < 1
    error('AK_GABAinASD_makeEyetrackingTableFigures requires you to designate a spreadsheet by name (string)')
end
if nargin < 2
    highlight_subjects = [];
    xlsfile_dir = 'L:\MurrayLab\ASD\SuppressionSummation\fMRI_HeadMotion_Data_Tables.xlsx';
end
if nargin < 3
    xlsfile_dir = 'L:\MurrayLab\ASD\SuppressionSummation\fMRI_HeadMotion_Data_Tables.xlsx';
end

% make sure highlight_subjects are cell arrays of stings
if ~isempty(highlight_subjects) && ~iscell(highlight_subjects)
    highlight_subjects = {highlight_subjects};
end

%% load eye tracking data tables

% read in data
[~,~,dataTable] = xlsread(xlsfile_dir,sheetName); % read contrast sheet

%% get dataTable column names and indices

% remove columns that don't correspond to data
removeColumns = ~cellfun(@ischar,dataTable(1,:));
dataTable(:,removeColumns) = [];

% find column indices for each stat to plot
for iPS = 1:length(dataTable(1,3:end)); 
    stat(iPS).name = dataTable{1,iPS+2};
    stat(iPS).Cindex = find(~cellfun(@isempty,regexp(dataTable(1,:),regexptranslate('wildcard',['*' dataTable{1,iPS+2}])))==1);
end

%% separate into vectors of data by group and stats of interest, generate labels, and group information into plotting function inputs

% fill in empty cells with nans
dataTable(cellfun(@isempty,dataTable)) = {nan};
dataTable = dataTable(2:end,:); % remove header row with labels

% % indices for group membership
% index for subjects of interest (highlight_subjects)
interestIndex = zeros(length(dataTable(:,1)),1); % preallocate
for iSubj = 1:length(highlight_subjects);
    interestIndex = interestIndex+strcmp(dataTable(:,1),highlight_subjects{iSubj});
end

% assignment of data arrays
for iStat = 1:length(stat)
    stat(iStat).values = arrayfun(@cell2mat,dataTable(:,stat(iStat).Cindex)); % for ASD group
end

% create labels for subject,fMRI session, and set#
 labels = cell(length(dataTable(:,1)),1); % preallocate for speed
 emptyLabel = {[]};
for iT = 1:length(dataTable(:,1))
    labels{iT} = [dataTable{iT,1}(3:end) ': set' num2str(dataTable{iT,2})];
end

% x-axis labels
x_labels = {'All Participants'};
 
% group data arrays, datapoint labels, and interest indices for plotting
if ~isempty(highlight_subjects)
    interestIndices = [];
end
for iStat = 1:length(stat)
    for iX = 1:length(x_labels);
        % match NT x_labels for cond with proper data arrays and empty datapoint labels
        stat(iStat).plotData(iX) = {stat(iStat).values};
        plotLabels{iX} = labels;
        % append interestIndices
        if ~isempty(highlight_subjects) && iStat==1
           interestIndices = vertcat(interestIndices,interestIndex);
        end
    end
end

%% figures

if ~isempty(highlight_subjects)
    for iStat = 1:length(stat)
        figure
        stat(iStat).handle = AK_plotSpread_mouseover(stat(iStat).plotData,plotLabels,'categoryIdx',interestIndices,'categoryColors',{'b','g'},'xNames',labels,'showMM',5);
        title(stat(iStat).name)
    end
else
    for iStat = 1:length(stat)
        figure
        stat(iStat).handle = AK_plotSpread_mouseover(stat(iStat).plotData,plotLabels,'xNames',labels,'showMM',5);
        title(stat(iStat).name)
    end
end

%% prepare output struct

for iStat = 1:length(stat)
    structure(iStat).statName = stat(iStat).name;
    structure(iStat).handle = stat(iStat).handle;
    structure(iStat).plotData = stat(iStat).plotData;
    structure(iStat).x_labels = x_labels;
    structure(iStat).mouseover_labels = plotLabels;
end
end

