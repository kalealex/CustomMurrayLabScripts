function [ FD, FDrms, FDsd, FDmax, setlist ] = AK_GABAinASD_getFDinfMRI( subject, fMRIsession, root_dir )
%AK_GABAinASD_getFDinfMRI accepts a subject code and fMRI session from the
%GABAinASD study as input and returns summary statistics on head motion for
%each set of functional data. The summary statistics returned focus on
%framewise displacement (FD) from Power et al. (2012).
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
%       FD: a cell array containing arrays of framewise displacements for
%           each set in the given session, where FD is the sum of the
%           absolute values of the displacements in each of six dimensions
%           of head motion (dx, dy, dz, rx, ry, rz) from one frame to the
%           next; rotional displacements converted from degrees to
%           millimeters by calculating displacements on the surface of a 
%           sphere with radius = 50 mm (Power et al. 2012); the length of
%           each set's array depends on the number of TRs in the scan
%       FDrms: an array containg the root mean squares of FDs for each set
%           in the given session
%       FDsd: an array containg the standard deviations of FDs for each set
%           in the given session
%       FDmax: an array containg the maximum values of FDs for each set in 
%           the given session
%       setlist: a cell array of strings naming the sets which correspond
%           to the other output arrays
%
% Reference
% Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage, 59(3), 2142–2154. http://doi.org/10.1016/j.neuroimage.2011.10.018


% check inputs
if nargin < 2
    error('AK_GABAinASD_getMotionCorrection requires at least two inputs: subject code and fMRI session name (both as strings)') % error message
end
if nargin < 3
    root_dir = 'L:\MurrayLab\ASD\Data';
end

% create home_dir
home_dir = fullfile(root_dir,subject,fMRIsession);
chdir(home_dir)

set_dirs = {'set1';'set2';'set3';'set4';'set5';'set6';'set7';'set8';'set9'}; % sets to check

%% Get motion parameters and Convert to FD (frame displacement): Based on Power et al 2012

% preallocate for speed
MC = cell(length(set_dirs),1);
FD = cell(length(set_dirs),1);
FDrms = nan(length(set_dirs),1);
FDsd = nan(length(set_dirs),1);
FDmax = nan(length(set_dirs),1);
setlist = [];

for iS = 1:length(set_dirs)
    if exist(fullfile(home_dir,set_dirs{iS}),'dir') % does directory exist?
        chdir(set_dirs{iS})
        if exist(sprintf('%s_3DMC.sdm', set_dirs{iS}),'file') % does file exist?
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
            FD{iS} = sum(FD{iS},2);

            % RMS of FD across scan to get 1 number to characterize scan motion
            FDrms(iS) = rms(FD{iS}); 
            % SD of FD
            FDsd(iS) = std(FD{iS}); 
            % maximum value of FD
            FDmax(iS) = max(FD{iS}); 

            % create list of sets included
            setlist = [setlist;set_dirs(iS)];
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

