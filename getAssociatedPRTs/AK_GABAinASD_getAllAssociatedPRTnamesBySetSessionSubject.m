function [ prtList ] = AK_GABAinASD_getAllAssociatedPRTnamesBySetSessionSubject( subjects, newListYN, listFilename )
%AK_GABAinASD_getAllAssociatedPRTnamesBySetSessionSubject puts together a
%list of subjects, fMRI sessions, sets and the prt files associated with them.
%This is useful for determining the correspondence between a given set and
%the type of scan that was run (i.e., contrast) and thus, the logfiles and
%condtion names writen to the eye tracking data files.
%   INPUT:
%       subjects: a cell array of strings containing a list of subject
%           codes to query; these may already be in an existing prtList
%       newListYN: a boolean where true tells the function to create a new
%           list of associated prt files and false tells the function to
%           compare your subjects list to the current copy of the list file
%           and append to that list only subjects who are not already there
%       listFilename: the name of the mat file containing the prtList;
%           defaults to 'AllAssociatedPRTfilenames.mat' which lives on the L
%           drive in 'L:\MurrayLab\ASD\Data'
%   OUTPUT:
%       prtList: a cell array containing the columns: subject, session,
%           set#, prtFilename; it has number of row equal to the number of
%           sets with existing directories
%
%   NOTES FOR USE: 
%       This function will only find associated prt file names for sets
%       where the functional data has been analyzed and linked to a prt
%       file. This means that unanalyzed data will not appear in prtList;
%       it also means that in order to add previously unanalyzed data to
%       the list for a subject who alread has at least one row in the list,
%       you must create a completely new list by changing the second input
%       to true. When you create a new list SAVE A COPY OF THE OLD LIST
%       because THE FUNCTION WILL NOT BE ABLE TO FIND ASSOCIATED PRT FILE
%       NAMES FOR SETS WITH BROKEN LINKS to prt files. This means that for 
%       these few sets, the names of associated prt files will have to be
%       manually entered. Saving a copy of the previous list makes this
%       manual entry as easy as copying and pasting into the blank cells in
%       the fourth column of the new prtList.

% directory info
top_dir = 'L:\MurrayLab\ASD\Data'; % set up base directory
fMRI_dirs = {'fMRI1','fMRI2','ftap'}; % designate folders for fMRI sessions
scan_dirs = {'set1','set2','set3','set4','set5','set6','set7','set8','set9'}; % designate folders for scans

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
    prtList = cell(length(subjects)*length(fMRI_dirs)*length(scan_dirs),4); % preallocate size
    prtList(1,:) = {'subject','session','set#','prtFilename'};
    
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
        prtList = cell(length(subjects)*length(fMRI_dirs)*length(scan_dirs),4); % preallocate size
        prtList(1,:) = {'subject','session','set#','prtFilename'};

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
            for iF = 1:length(fMRI_dirs) % cycle through sessions
                for iScan = 1:length(scan_dirs) %cycle through scans
                    use_dir = fullfile(subj_dir,fMRI_dirs{iF},scan_dirs{iScan}); % generate directory in use
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
                        prtList(prtListRow,:) = {subj_dirs{iS},fMRI_dirs{iF},str2double(scan_dirs{iScan}(4)),PRTfilename};
                        prtListRow = prtListRow+1; % advance counter


                    end
                end
            end
        end
    end
    % trim empty rows from cell array
    prtList(all(cellfun(@isempty,prtList(:,:)),2),:) = [];
end

save(listFilename,'prtList'); % save the current list

end

