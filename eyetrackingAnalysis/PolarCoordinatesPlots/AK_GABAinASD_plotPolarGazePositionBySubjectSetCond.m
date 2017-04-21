function [ h ] = AK_GABAinASD_plotPolarGazePositionBySubjectSetCond( subjects, sheets, dcCoordinates, nttCriterion, xlsfile_dir )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%% inputs

% temp; to be selected in GUI
sheets = {'Contrast','V1_fix Loc'}; % which sheets
subjects = {'G114','G312'}; % which subjects
dcCoordinates = 1; % drift corrected coordinates?
nttCriterion = 1; % filter out data with no tracking time values larger than this criterion
xlsfile_dir = 'L:\MurrayLab\ASD\Data\fMRI_EyeTracking_Data_Tables_OnlyGood.xlsx'; % file to load eye tracking data from 

% drift corrected or raw?
if dcCoordinates==1
    plotStats = {'drift corrected distance M','drift corrected angle M','drift corrected distance SD','drift corrected angle SD'};
else
    plotStats = {'distance M','angle M','distance SD','angle SD'};
end

%% load eye tracking data tables for each sheetl; extract polar coordinates and labels by condition
                
% preallocate
subj(1:length(subjects)) = struct('dist',[],'angle',[],'condIdx',[],'mLabelIdx',[],'legendLabels',[],'mouseoverLabels',[]); % struct for plot data
condsPerSheet = zeros(length(sheets),1); % for counting number of conditions across sheets
setsPerSS = zeros(length(sheets),length(subjects)); % for counting number of sets across sheets and subjects
plotDataIdx = ones(length(subjects),1); % set counter for plotData arrays

for iSheet = 1:length(sheets)
% load sheet and determine whether it is a block sheet:    
    
    % this sheet
    clear sheetName dataTable isBlockSheet nttColumnIdx stat findCondStr condStrIndex fullCondList condList condIdx
    sheetName = sheets{iSheet};
    
    % read in data
    [~,~,dataTable] = xlsread(xlsfile_dir,sheetName);

    % is this a sheet for block statistics?
    isBlockSheet = ~isempty(regexp(sheetName,regexptranslate('wildcard','block*'),'ONCE'));

% create dataTable column indices from list of plotStats; get condition names:

    % find column indices for no tracking time
    nttColumnIdx = find(~cellfun(@isempty,regexp(dataTable(1,:),regexptranslate('wildcard','*no tracking time')))==1);

    % find column indices for each stat to plot
    for iPS = 1:length(plotStats); 
        stat(iPS).name = plotStats{iPS};
        stat(iPS).Cindex = find(~cellfun(@isempty,regexp(dataTable(1,:),regexptranslate('wildcard',['*' plotStats{iPS}])))==1);
    end

    % find condition names
    findCondStr = '*:*';
    condStrIndex = find(~cellfun(@isempty,regexp(dataTable(1,:),regexptranslate('wildcard',findCondStr)))==1);
    for iCS = 1:length(condStrIndex)
        clear tempIdx
        tempIdx = strfind(dataTable{1,condStrIndex(iCS)},findCondStr(2)); % idx for position of : in each string where : occurs
        fullCondList{iCS} = dataTable{1,condStrIndex(iCS)}(1:tempIdx-1); % list of all condition names in order of columns
    end
    condList = unique(fullCondList,'stable');
    condsPerSheet(iSheet) = length(condList); % count conds per sheet
    % creat column index for each condition
    for iCL = 1:length(condList);
        condIdx{iCL} = find(strcmp(fullCondList,condList{iCL})) + condStrIndex(1) - 1; % create an index for columns for each condition
    end
    
% filter by nttCriterion; create row indices for subjects; separate data into vectors of data by subject, stats and cond; generate legend labels and mouseover labels:

    % this step is different for block and summary stats
    switch isBlockSheet
        case 0 % for summary stats:
            clear nttFilterIndex subjIdx
            
            % fill in empty cells with nans
            dataTable(cellfun(@isempty,dataTable)) = {nan};

            % filter out values for scans with no track time > criterion 
            nttFilterIndex = [nan;arrayfun(@cell2mat,dataTable(2:end,nttColumnIdx))] < nttCriterion; % these are rows to keep
            dataTable = dataTable(nttFilterIndex,:);
            
            for iSubj = 1:length(subjects)
                % find row indices for subjects
                subjIdx{iSubj} = find(strcmp(dataTable(:,1),subjects{iSubj}));
                
                % document #sets in this sheet for this subject
                setsPerSS(iSheet,iSubj) = length(subjIdx{iSubj});
                
                for iSet = 1:length(subjIdx{iSubj})
                    % assignment of data arrays
                    for iCond = 1:length(condList)
                        % joint index for stat and cond
                        clear distMCondIdx angleMCondIdx distSDCondIdx angleSDCondIdx
                        distMCondIdx = intersect(stat(1).Cindex,condIdx{iCond}); % should be only one column where stat and cond intersect
                        angleMCondIdx = intersect(stat(2).Cindex,condIdx{iCond});
                        distSDCondIdx = intersect(stat(3).Cindex,condIdx{iCond}); 
                        angleSDCondIdx = intersect(stat(4).Cindex,condIdx{iCond});
                        % store column values for dist and angle by subject
                        subj(iSubj).distM(plotDataIdx(iSubj)) = cell2mat(dataTable(subjIdx{iSubj}(iSet),distMCondIdx));
                        subj(iSubj).angleM(plotDataIdx(iSubj)) = cell2mat(dataTable(subjIdx{iSubj}(iSet),angleMCondIdx));
                        subj(iSubj).distSD(plotDataIdx(iSubj)) = cell2mat(dataTable(subjIdx{iSubj}(iSet),distSDCondIdx));
                        subj(iSubj).angleSD(plotDataIdx(iSubj)) = cell2mat(dataTable(subjIdx{iSubj}(iSet),angleSDCondIdx));
                        % store condition and mouseoverLabel Indices
                        if iSheet==1 
                            subj(iSubj).condIdx(plotDataIdx(iSubj)) = iCond;
                            subj(iSubj).mLabelIdx(plotDataIdx(iSubj)) = plotDataIdx(iSubj);
                        else
                            subj(iSubj).condIdx(plotDataIdx(iSubj)) = sum(condsPerSheet(1:iSheet-1)) + iCond;
                            subj(iSubj).mLabelIdx(plotDataIdx(iSubj)) = plotDataIdx(iSubj);
                        end
                        plotDataIdx(iSubj) = plotDataIdx(iSubj)+1; % advance array counter
                        
                        % create labels for legend
                        if iSet==1 % avoid doing this more than once per condition and sheet
                            subj(iSubj).legendLabels = [subj(iSubj).legendLabels {[sheetName ': ' condList{iCond}]}];
                        end
                        % create mouseover labels for subject,fMRI session, and set#
                        subj(iSubj).mouseoverLabels = [ subj(iSubj).mouseoverLabels {sprintf([dataTable{subjIdx{iSubj}(iSet),2} ': set' num2str(dataTable{subjIdx{iSubj}(iSet),3}) '\n'...
                            'distSD: ' num2str(cell2mat(dataTable(subjIdx{iSubj}(iSet),distSDCondIdx))) '\n'...
                            'angleSD: ' num2str(cell2mat(dataTable(subjIdx{iSubj}(iSet),angleSDCondIdx)))])} ];
                        
                    end
                    
                end    
            end

        case 1 % for block stats:
%             error('This function is not yet developed for plotting block stats') % message

            % fill in empty cells with nans
            dataTable(cellfun(@isempty,dataTable)) = {nan};

            % filter out values for scans with no track time > criterion 
            for iCnd = 1:length(condList);
                clear nttFilterIndex
                nntFilterIndex = [nan;arrayfun(@cell2mat,dataTable(2:end,nttColumnIdx(iCnd)))] >= nttCriterion; % these are rows to remove
                dataTable(nntFilterIndex,condIdx{iCnd}) = {nan}; % nans and corresponding labels will be removed by AK_compass_mouseover
            end
            dataTable = dataTable(2:end,:); % remove header row with labels
        
            for iSubj = 1:length(subjects)
                % find row indices for subjects
                subjIdx{iSubj} = find(strcmp(dataTable(:,1),subjects{iSubj}));
                
                % document #sets in this sheet for this subject (really number of blocks)
                setsPerSS(iSheet,iSubj) = length(subjIdx{iSubj});
                
                for iSet = 1:length(subjIdx{iSubj}) % this cycles through each block
                    % assignment of data arrays
                    for iCond = 1:length(condList)
                        % joint index for stat and cond
                        clear distMCondIdx angleMCondIdx distSDCondIdx angleSDCondIdx
                        distMCondIdx = intersect(stat(1).Cindex,condIdx{iCond}); % should be only one column where stat and cond intersect
                        angleMCondIdx = intersect(stat(2).Cindex,condIdx{iCond});
                        distSDCondIdx = intersect(stat(3).Cindex,condIdx{iCond}); 
                        angleSDCondIdx = intersect(stat(4).Cindex,condIdx{iCond});
                        % store column values for dist and angle by subject
                        subj(iSubj).distM(plotDataIdx(iSubj)) = cell2mat(dataTable(subjIdx{iSubj}(iSet),distMCondIdx));
                        subj(iSubj).angleM(plotDataIdx(iSubj)) = cell2mat(dataTable(subjIdx{iSubj}(iSet),angleMCondIdx));
                        subj(iSubj).distSD(plotDataIdx(iSubj)) = cell2mat(dataTable(subjIdx{iSubj}(iSet),distSDCondIdx));
                        subj(iSubj).angleSD(plotDataIdx(iSubj)) = cell2mat(dataTable(subjIdx{iSubj}(iSet),angleSDCondIdx));
                        % store condition and mouseoverLabel Indices
                        if iSheet==1 
                            subj(iSubj).condIdx(plotDataIdx(iSubj)) = iCond;
                            subj(iSubj).mLabelIdx(plotDataIdx(iSubj)) = plotDataIdx(iSubj);
                        else
                            subj(iSubj).condIdx(plotDataIdx(iSubj)) = sum(condsPerSheet(1:iSheet-1)) + iCond;
                            subj(iSubj).mLabelIdx(plotDataIdx(iSubj)) = plotDataIdx(iSubj);
                        end
                        plotDataIdx(iSubj) = plotDataIdx(iSubj)+1; % advance array counter
                        
                        % create labels for legend
                        if iSet==1 % avoid doing this more than once per condition and sheet
                            subj(iSubj).legendLabels = [subj(iSubj).legendLabels {[sheetName ': ' condList{iCond}]}];
                        end
                        
                        % create mouseover labels for subject,fMRI session, and set#
                        subj(iSubj).mouseoverLabels = [ subj(iSubj).mouseoverLabels {sprintf([dataTable{subjIdx{iSubj}(iSet),2} ': set' num2str(dataTable{subjIdx{iSubj}(iSet),3}) ': block' num2str(dataTable{subjIdx{iSubj}(iSet),5}) '\n'...
                            'distSD: ' num2str(cell2mat(dataTable(subjIdx{iSubj}(iSet),distSDCondIdx))) '\n'...
                            'angleSD: ' num2str(cell2mat(dataTable(subjIdx{iSubj}(iSet),angleSDCondIdx)))])} ];
                    end
                    
                end    
            end
    end
end

% generate colors for table
condCount = sum(condsPerSheet);
plotColors = distinguishable_colors(condCount);

%% figures

for iSu = 1:length(subjects)
    % new figure and axis for each subject
    h(iSu).Fig = figure; 
    h(iSu).Ax = axes('Parent',h(iSu).Fig,'visible','on');
    
    % plot colors offscreen for legend
    for iColor = 1:length(plotColors(:,1))
        line([1000 1000],[1000+iColor 1000],'Color',plotColors(iColor,:)); hold on;
    end
    
    % plot patches to represent SD in angle and distance
    h(iSu).Patch = AK_polarPatch(h(iSu).Ax,subj(iSu).distM,subj(iSu).angleM,subj(iSu).distSD,subj(iSu).angleSD,plotColors,subj(iSu).condIdx);
    
    % plot and store handles
    h(iSu).Plot = AK_compass_mouseover(h(iSu).Ax,subj(iSu).distM,subj(iSu).angleM,plotColors,subj(iSu).condIdx,subj(iSu).mouseoverLabels,subj(iSu).mLabelIdx);
    
    % legend
    h(iSu).Leg = legend(subj(iSu).legendLabels);
      
    %title
    clear TitleStr
    TitleStr = subjects{iSu};
    for iSh = 1:length(sheets)
        if iSheet==1
            TitleStr = [TitleStr ': ' sheets{iSh}];
        else
            TitleStr = [TitleStr ', ' sheets{iSh}];
    end
    title(TitleStr);
end

end

