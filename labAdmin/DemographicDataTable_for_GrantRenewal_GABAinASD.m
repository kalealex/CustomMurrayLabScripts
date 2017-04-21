% DemographicDataTable_for_GrantRenewal_GABAinASD
% Murray Lab 2016
% created by AMK 12/19/16

%% directories and demographic data file name

root_dir = 'C:\Users\Alex Kale\Documents\MATLAB\MurrayLab\GrantRenewalDemographicData';

demo_filename = fullfile(root_dir,'DemographicData161219_ExcludeDQ.xlsx');

%% format table for Scott

rowLabels = {'Native American';'Asian';'Pacific Islander';'Black';'White';'Multiracial';'Unknown'};
columnLabels = {'Not Hispanic: Female','Not Hispanic: Male','Not Hispanic: Unknown','Hispanic: Female','Hispanic: Male','Hispanic: Unknown','Not Reported: Female','Not Reported: Male','Not Reported: Unknown'};
demograpicsMatrix = zeros(length(rowLabels),length(columnLabels));

%% load file containing demographic info

% load
[~,~,demoCell] = xlsread(demo_filename);

% setup
NoExp = regexptranslate('wildcard','Unchecked*'); % for race questions
YesExp = regexptranslate('wildcard','Checked*');
nConsentedButNotIncluded = 0; % counter

% fill in demographics matrix from demographics data in .xls file
for iD = 3:length(demoCell(:,1)) % start at third row
    % count participants not included in matrix 
    % no data entered for subjects with incomplete data on REDCap
    if isnan(demoCell{iD,2}) && logical(mean(isnan(demoCell{iD,10}))) && all(~cellfun(@isempty,cellfun(@(x) regexp(x,NoExp),demoCell(iD,3:8),'UniformOutput',false))) 
       nConsentedButNotIncluded = nConsentedButNotIncluded+1;
    else
        % figure out which column and row to add this participant to
        clear addCol addRow
        % columns...
        if ~isempty(regexp(demoCell(iD,10),regexptranslate('wildcard','No*'),'once')) % is not hispanic? (strcmp wasn't working here)
            if strcmp(demoCell(iD,2),'F') % is female?
                addCol = 1;
            elseif strcmp(demoCell(iD,2),'M') % is male?
                addCol = 2;
            elseif isnan(demoCell{iD,2}) % is unknown?
                addCol = 3;
                disp(['subject ' demoCell{iD,1} ' no gender']);
            end
        elseif ~isempty(regexp(demoCell(iD,10),regexptranslate('wildcard','Yes*'),'once')) % is hispanic? (strcmp wasn't working here)
            if strcmp(demoCell(iD,2),'F') % is female?
                addCol = 4;
            elseif strcmp(demoCell(iD,2),'M') % is male?
                addCol = 5;
            elseif isnan(demoCell{iD,2}) % is unknown?
                addCol = 6;
            end
        elseif logical(mean(isnan(demoCell{iD,10}))) % is not reported hispanic?
            if strcmp(demoCell(iD,2),'F') % is female?
                addCol = 7;
            elseif strcmp(demoCell(iD,2),'M') % is male?
                addCol = 8;
            elseif isnan(demoCell{iD,2}) % is unknown?
                addCol = 9;
            end    
        end
        % rows... (strcmp wasn't working here)
        if sum(~cellfun(@isempty,cellfun(@(x) regexp(x,YesExp),demoCell(iD,3:8),'UniformOutput',false))) > 1 || strcmp(demoCell(iD,9),'Multi-ethnic') % is multiracial?
            addRow = 6;
        elseif all(~cellfun(@isempty,cellfun(@(x) regexp(x,NoExp),demoCell(iD,3:8),'UniformOutput',false))) % is unknown?
            addRow = 7;
        elseif ~cellfun(@isempty,regexp(demoCell(iD,3),YesExp,'once')) % is asian?
            addRow = 2;
        elseif ~cellfun(@isempty,regexp(demoCell(iD,4),YesExp,'once')) % is black?
            addRow = 4; 
        elseif ~cellfun(@isempty,regexp(demoCell(iD,5),YesExp,'once')) % is native american?
            addRow = 1;  
        elseif ~cellfun(@isempty,regexp(demoCell(iD,6),YesExp,'once')) % is islander?
            addRow = 3;
        elseif ~cellfun(@isempty,regexp(demoCell(iD,7),YesExp,'once')) || strcmp(demoCell(iD,9),'Middle Eastern') % is white?
            addRow = 5;    
        end
        % add subject to correct position in matrix
        demograpicsMatrix(addRow,addCol) = demograpicsMatrix(addRow,addCol) + 1;
    end  
end

%% save

save(fullfile(root_dir,'DemographicData_forScott.mat'),'demograpicsMatrix','rowLabels','columnLabels','nConsentedButNotIncluded')
