function [ motionCorrection ] = AK_read3DMC( file_dir )
%AK_read3DMC reads BVQX 3D motion correction output files named *3DMC.log
%and returns slice-by-slice motion correction data in six dimensions 
%   INPUT:
%       file_dir: a full file name (string) for a formatted *3DMC.log file;
%           these files are created by the BVQX motion correction algorithm
%           we use for GABAinASD.
%   OUTPUT:
%       motionCorrection: a n-by-m matrix where n is equal to the number of
%           TRs and m represents the following six dimensions in order: 
%           dx (mm), dy (mm), dz (mm), rx (deg), ry (deg), rz (deg) 

% check input
if nargin ~= 1 || ~ischar(file_dir) || isempty(regexp(file_dir,regexptranslate('wildcard','*3DMC.log'),'once'))
    error(['AK_read3DMC requires a full file name (string) for a formatted *3DMC.log file as input. '
        'These files are created by the BVQX motion correction algorithm we use for GABAinASD.']);
end

% read file
rawText = fileread(file_dir);
% split into rows
textRowCell = strsplit(rawText,'-> ');
textRowCell(1) = []; % remove header
% parse each row into values of interest
motionCorrection = nan(length(textRowCell),6); % preallocate
for iTR = 1:length(textRowCell) % each row represents a TR
    motionCorrection(iTR,1) = getValueAfter(textRowCell{iTR},'dx: ');
    motionCorrection(iTR,2) = getValueAfter(textRowCell{iTR},'dy: ');
    motionCorrection(iTR,3) = getValueAfter(textRowCell{iTR},'dz: ');
    motionCorrection(iTR,4) = getValueAfter(textRowCell{iTR},'rx: ');
    motionCorrection(iTR,5) = getValueAfter(textRowCell{iTR},'ry: ');
    motionCorrection(iTR,6) = getValueAfter(textRowCell{iTR},'rz: ');
end

    % gets the full numeric value from the string myStr which immediately follows the substring mySubstr
    function [value] = getValueAfter(myStr, mySubstr)
        mySubstrIdx = strfind(myStr,mySubstr) + length(mySubstr); % starting index where number is expected
        numCharCount = 0; 
        while (~isspace(myStr(mySubstrIdx + numCharCount)))
            numCharCount = numCharCount + 1;
        end
        value = str2double(myStr(mySubstrIdx:mySubstrIdx + numCharCount - 1));
    end

end

