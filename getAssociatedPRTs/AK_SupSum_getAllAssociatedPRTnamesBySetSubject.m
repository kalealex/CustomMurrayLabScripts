function [ prtList ] = AK_SupSum_getAllAssociatedPRTnamesBySetSubject( subjects, newListYN, listFilename )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% directory info
top_dir = 'L:\MurrayLab\ASD\SuppressionSummation'; % set up base directory
scan_dirs = {'set1','set2','set3','set4','set5','set6','set7','set8'}; % designate folders for scans

% check inputs
if nargin < 1 || ~iscell(subjects)
    error('AK_GABAinASD_getAllAssociatedPRTnamesBySetSessionSubject requires a list of subject IDs as a cell array of strings')
elseif nargin < 2 
    newListYN = 0; % default to appending onto existing list
elseif nargin < 3
    listFilename = fullfile(top_dir,'AllAssociatedPRTfilenames.mat');
end

if newListYN==1 % prepare a new list
    % prepare output cell array
    prtList = cell(length(subjects)*length(scan_dirs),3); % preallocate size
    prtList(1,:) = {'subject','set#','prtFilename'};
    
    subj_dirs = subjects; % include all subjects in subj_dirs
    
    prtListRow = 2; % create row counter
else
    if exist(listFilename,'file') % check that file exists
        load(listFilename,'prtList'); % load existing list

        % find list of subjects not already included in the cell array
        prtListSubjects = unique(prtList(2:end,1)); % list of unique subject names in existing list
        subj_dirs = setdiff(subjects,prtListSubjects); % make subj_dirs = list of previously unincluded subjects

        prtListRow = length(prtList(:,1))+1; % create row counter  
    else
        disp('***Could not find prtList file; creating new list') % message
        % prepare output cell array
        prtList = cell(length(subjects)*length(scan_dirs),3); % preallocate size
        prtList(1,:) = {'subject','set#','prtFilename'};

        subj_dirs = subjects; % include all subjects in subj_dirs

        prtListRow = 2; % create row counter
    end
end

if ~isempty(subj_dirs) % check to make sure there are subjects who need to be added to the list
    for iS = 1:length(subj_dirs) % cycle through subjects
        subj_dir = fullfile(top_dir,subj_dirs{iS});
        if exist(subj_dir,'dir')==0
            disp(['could not find subject folder for' subj_dirs{iS}]) % message
        elseif exist(subj_dir,'dir')==7
            for iScan = 1:length(scan_dirs) %cycle through scans
                use_dir = fullfile(subj_dir,scan_dirs{iScan}); % generate directory in use
                if exist(use_dir,'dir')==0 
                    disp([use_dir ' is not an existing directory']) % message
                else
                    % create setIDfilename name to look for previously saved prt file name
                    clear setIDfilename PRTfilename
                    setIDfilename = [fullfile(use_dir,scan_dirs{iScan}) '_scanID.mat'];

                    % either load saved filename from mat file or query filename using BVQXfile
                    try 
                        if exist(setIDfilename,'file')==0; % does the set scanID file exist?
                            clear vtcfileName vtc 
                            BVQXfile(0, 'clearallobjects');
                            vtcfileName = [fullfile(use_dir,scan_dirs{iScan}) '.vtc']; % name vtc files 
                            vtc = BVQXfile(deblank(vtcfileName)); % fill out vtc structure for purpose of naming scans
                            PRTfilename = vtc.NameOfLinkedPRT;
                            save(setIDfilename,'PRTfilename'); % save associated .prt filename
                        else
                            load(setIDfilename,'PRTfilename'); % load associated .prt filename
                        end
                    catch
                        PRTfilename = [];
                        disp(['trouble processing ' vtcfileName]); % message
                    end

                    % append filename by set to cell array
                    prtList(prtListRow,:) = {subj_dirs{iS},str2double(scan_dirs{iScan}(4)),PRTfilename};
                    prtListRow = prtListRow+1; % advance counter


                end
            end
        end
    end
    % trim empty rows from cell array
    prtList(all(cellfun(@isempty,prtList(:,:)),2),:) = [];
end

save(listFilename,'prtList'); % save the current list

end

