function [ set, FDrms, FDsd, FDmax ] = AK_GABAinASD_fMRI_getFramewiseDisplacementPerCondition( subject, fMRIsession, root_dir )
%AK_GABAinASD_fMRI_getFramewiseDisplacementPerCondition accepts a subject
%code and fMRI session from the GABAinASD study as input and returns
%summary statistics on head motion for each set of functional data. The
%summary statistics returned focus on framewise displacement (FD) from
%Power et al. (2012). The function also returns a structure that contains
%framewise displacements for each TR, indices identifying the stimulus
%condition for each TR, the experiment, and the directory for each set of
%fMRI.
%   INPUT: 
%       subject: a subject code from the GABAinASD study (string)
%       fMRIsession: an fMRI sessionn name from the GABAinASD study (string); 
%           either 'fMRI1', 'fMRI2' or 'ftap' (a list of sets checked
%           against the existing directory structure)
%       root_dir[optional]: the root directory containing the subject
%           folders from the GABAinASD study; should always be the
%           default value 'L:\MurrayLab\ASD\Data', so no need for this
%           input unless you want to look at a local copy of the data
%   OUTPUT:
%       set: a structure with one unit of length per available set of fMRI
%           data in the given fMRI session; the structure contains the
%           following fields:
%           directory: the directory where the data about this set can be
%               found
%           FD: an array of framewise displacements (FD), where FD is the
%               sum of the absolute values of the displacements in each of
%               six dimensions of head motion (dx, dy, dz, rx, ry, rz) from
%               one frame (TR) to the next; rotional displacements
%               converted from degrees to millimeters by calculating
%               displacements on the surface of a sphere with radius = 50
%               mm (Power et al. 2012); the length the array is equal to
%               the number of TRs in the scan, and the first FD in the list
%               will always be zero
%           experiment: a string identifying the fMRI experiment which was
%               run during a given set; selected from the list:
%               'MTlocalizer' ,'contrast', 'suppression', 'summation',
%               'V1localizer', 'V1localizer_fix', 'ftap1', 'ftap2', 'ftap3'
%           conditionIdx: an index into the array FD where each index value
%               represents a different condition from the fMRI experiment
%               (see the field conditionList, which is a key for this
%               index)
%           conditionList: a cell array of strings, a key for the
%               conditionIdx field; the values of the indices in
%               conditionIdx match the positions in this cell array of the
%               strings naming each condition in this experiment 
%       FDrms: an array containg the root mean squares of FDs for each set
%           in the given session
%       FDsd: an array containg the standard deviations of FDs for each set
%           in the given session
%       FDmax: an array containg the maximum values of FDs for each set in 
%           the given session
%
% Reference
% Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage, 59(3), 2142–2154. http://doi.org/10.1016/j.neuroimage.2011.10.018
%
% History: This is an updated version of the function
% AK_GABAinASD_getFDinfMRI with the added feature that the function queries
% a stored list of linked .prt file names for each set in the given session
% in order to provide a link between individual TRs in each set and the
% stimulus conditions to which they correspond. The purpose of this change
% is to enable the flexible computation of summary statistics on FD at the
% level of conditions per run, or even at the level of individual blocks of
% TRs, such that it is possible to identify fMRI data which may be
% corrupted by head motion and remove it from further analysis.



% check inputs
if nargin < 2
    error('AK_GABAinASD_getMotionCorrection requires at least two inputs: subject code and fMRI session name (both as strings)') % error message
end
if nargin < 3
    root_dir = 'L:\MurrayLab\ASD\Data';
end


%% get a list of associated .prt file names for this subject (used to parse head motion by conditions)

% load existing .prt file list
prtList = AK_GABAinASD_getAllAssociatedPRTnamesBySetSessionSubject({subject},false);
% reduce list to requisite information
prtList = prtList(strcmp(prtList(:,1),subject) & strcmp(prtList(:,2),fMRIsession),3:4);
% match formatting of set numbers to set_dirs (see below)
prtList(:,1) = cellfun(@(x) ['set' num2str(x)],prtList(:,1),'UniformOutput',false);

%% hard code list of TR indices per condition

% cell array of possible .prt file names to match
experimentIDs = {'MTlocalizer','contrast','suppression','summation','V1localizer','V1localizer_fix','ftap1','ftap2','ftap3'}; % indices match table below

% table of condition names and TR indices (from .prt files) listed in order of experimentIDs
condIdxTable = {table([1:5, 11:15, 21:25, 31:35, 41:45, 51:55, 61:65, 71:75, 81:85, 91:95, 101:105, 111:115, 121:125],... % MT localizer
                    [6:10, 16:20, 26:30, 36:40, 46:50, 56:60, 66:70, 76:80, 86:90, 96:100, 106:110, 116:120],...
                    'VariableNames',{'static','moving'}),...
                table([1:5, 11:15, 21:25, 31:35, 41:45, 51:55, 61:65, 71:75, 81:85, 91:95, 101:105, 111:115, 121:125],... % contrast
                    [6:10, 26:30, 46:50, 66:70, 86:90, 106:110],...
                    [16:20, 36:40, 56:60, 76:80, 96:100, 116:120],...
                    'VariableNames',{'fixation','low','high'}),...
                table([1:5, 11:15, 21:25, 31:35, 41:45, 51:55, 61:65, 71:75, 81:85, 91:95, 101:105, 111:115, 121:125],... % suppression
                    [6:10, 16:20, 26:30, 36:40, 46:50, 56:60, 66:70, 76:80, 86:90, 96:100, 106:110, 116:120],...
                    'VariableNames',{'small','big'}),...
                table([1:5, 11:15, 21:25, 31:35, 41:45, 51:55, 61:65, 71:75, 81:85, 91:95, 101:105, 111:115, 121:125],... % summation (same as suppression)
                    [6:10, 16:20, 26:30, 36:40, 46:50, 56:60, 66:70, 76:80, 86:90, 96:100, 106:110, 116:120],...
                    'VariableNames',{'small','big'}),...
                table([1:5, 11:15, 21:25, 31:35, 41:45, 51:55, 61:65, 71:75],... % V1 localizer
                    [6:10, 16:20, 26:30, 36:40, 46:50, 56:60, 66:70, 76:80],...
                    'VariableNames',{'small','big'}),...
                table([1:5, 11:15, 21:25, 31:35, 41:45, 51:55, 61:65, 71:75],... % V1_fix localizer
                    [6:10, 16:20, 26:30, 36:40, 46:50, 56:60, 66:70, 76:80],...
                    'VariableNames',{'small','fixation'}),...
                table([1:10, 21:30, 41:50, 61:70, 81:90, 101:110, 121:130, 141:150, 161:170],... % ftap1
                    [11:20, 31:40, 111:120, 151:160],...
                    [51:60, 71:80, 91:100, 131:140],...
                    'VariableNames',{'fixation','standard','variable'}),...
                table([1:10, 21:30, 41:50, 61:70, 81:90, 101:110, 121:130, 141:150, 161:170],... % ftap2
                    [51:60, 111:120, 131:140, 151:160],...
                    [11:20, 31:40, 71:80, 91:100],...
                    'VariableNames',{'fixation','standard','variable'}),...
                table([1:10, 21:30, 41:50, 61:70, 81:90, 101:110, 121:130, 141:150, 161:170],... % ftap3
                    [11:20, 51:60, 71:80, 151:160],...
                    [31:40, 91:100, 111:120, 131:140],...
                    'VariableNames',{'fixation','standard','variable'})
                };
                

%% Get motion parameters and Convert to FD (frame displacement): Based on Power et al 2012

% create home_dir
home_dir = fullfile(root_dir,subject,fMRIsession);
chdir(home_dir)

% sets to check
set_dirs = {'set1';'set2';'set3';'set4';'set5';'set6';'set7';'set8';'set9'}; 

% preallocate for speed
MC = cell(length(set_dirs),1);
FD = cell(length(set_dirs),1);
set = struct;
FDrms = nan(length(set_dirs),1);
FDsd = nan(length(set_dirs),1);
FDmax = nan(length(set_dirs),1);

for iS = 1:length(set_dirs)
    clear thisSetDir
    thisSetDir = fullfile(home_dir,set_dirs{iS});
    if exist(thisSetDir,'dir') % does directory exist?
        chdir(set_dirs{iS})
        if exist(sprintf('%s_3DMC.log', set_dirs{iS}),'file') % does file exist?
            % read in motion correction coordinates
            evalStr = sprintf('[xt{iS} yt{iS} zt{iS} xr{iS} yr{iS} zr{iS}] = textread(''%s_3DMC.sdm'',''%s'',''headerlines'',9);', set_dirs{iS}, '%f %f %f %f %f %f');
            eval(evalStr)
            % reorganize into a single matrix array
            MC{iS} = [xt{iS} yt{iS} zt{iS} xr{iS} yr{iS} zr{iS}];

            % reset dir
            chdir(home_dir)

            % assign MC to tempory array
            clear temp
            temp = MC{iS}; 
            % convert degrees to radians to mm on a 50 mm circle
            temp(:,4:6) = ((temp(:,4:6).*pi)./180).*50;
            % difference between timepoint i and timepoint i-1
            FD{iS} = abs(temp(2:size(temp,1),:) - temp(1:size(temp,1)-1,:));
            % sum across all motion parameters
            FD{iS} = [0; sum(FD{iS},2)]; % add displacement of zero at the first index
            
            % select indices to separate framewise displacements by
            % condition based on associated .prt file name
            clear thisSetPRT thisSetExperimentIdx experimentConds condIdx
            thisSetPRT = prtList(strcmp(set_dirs{iS},prtList(:,1)),2);
            if ~isempty(thisSetPRT)
                thisSetExperimentIdx = ~cellfun(@isempty, cellfun(@(x) regexp(thisSetPRT,['.' x '\.prt']), experimentIDs));
                % create index of conditions at each TR
                experimentConds = condIdxTable{thisSetExperimentIdx}.Properties.VariableNames;
                condIdx = zeros(size(FD{iS})); % preallocate
                for iC = 1:length(experimentConds)
                    condIdx(condIdxTable{thisSetExperimentIdx}{1,iC}) = iC; % index numbers should match position in experimentConds; should match positions of FD{iS}
                end

                % organize FD data by condition into output structure
                set(iS).directory = thisSetDir;
                set(iS).experiment = experimentIDs{thisSetExperimentIdx};
                set(iS).FD = FD{iS};
                set(iS).conditionIdx = condIdx;
                set(iS).conditionList = experimentConds;
            else % no .prt file name for this set
                warning(['Missing .prt file link; cannot identify experiment for ' set_dirs{iS}]);
                
                % cannot organize FD data by condition; still create output structure
                set(iS).directory = thisSetDir;
                set(iS).experiment = 'cannot identify experiment because of missing .prt file link';
                set(iS).FD = FD{iS};
                set(iS).conditionIdx = nan;
                set(iS).conditionList = 'cannot identify condition set because of missing .prt file link';
            end
            
            % RMS of FD across scan to get 1 number to characterize scan motion
            FDrms(iS) = rms(FD{iS}); 
            % SD of FD
            FDsd(iS) = std(FD{iS}); 
            % maximum value of FD
            FDmax(iS) = max(FD{iS}); 
        else
            disp([sprintf('%s_3DMC.sdm', set_dirs{iS}) ' file does not exist for ' subject ': ' fMRIsession]) % message
        end
        % jump back home
        cd(home_dir);
    else
        disp([set_dirs{iS} ' directory does not exist for ' subject ': ' fMRIsession]) % message
    end
end


end


