function [ FD, FDrms, FDsd, FDmax, setlist ] = AK_SupSum_getFDinfMRI( subject, root_dir )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% check inputs
if nargin < 1
    error('AK_GABAinASD_getMotionCorrection requires at least one input: subject code (string)') % error message
end
if nargin < 2
    root_dir = 'L:\MurrayLab\ASD\SuppressionSummation';
end

% create home_dir
home_dir = fullfile(root_dir,subject);
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
    if exist(fullfile(home_dir,set_dirs{iS}),'dir')
        chdir(set_dirs{iS})
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
        disp([set_dirs{iS} ' does not exist for ' subject]) % message
    end
end



end

