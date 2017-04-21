function [ PUI, CM, PDRmin ] = AK_PupilSize_FatigueMetrics_Miosis( pupilSize, chunkSize )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% determine number of chunks; remainder of pupil signal will be unused
nChunks = idivide(length(pupilSize),int32(chunkSize));

% parse pupil size signal into nChunks and calculate miosis metrics:
for iChunk = 1:nChunks
    % find mean pupil size of chunk
    pupilM(iChunk) = nanmean(pupilSize(chunkSize*(iChunk-1)+1:chunkSize*iChunk));
    % find absolute differences between chunk means
    if iChunk > 1
        pupilMdiff(iChunk-1) = abs(pupilM(iChunk)-pupilM(iChunk-1));
    end
    % find ratio of each mean to the mean of the first chunk
    PDR = pupilM(iChunk)/pupilM(1);
end
% calculate PUI (pupillary unrest index)
    % Lüdtke, H., Wilhelm, B., Adler, M., Schaeffel, F., & Wilhelm, H. (1998). Mathematical procedures in data recording and processing of pupillary fatigue waves. Vision Research, 38(19), 2889–2896. http://doi.org/10.1016/S0042-6989(98)00081-9
PUI = sum(pupilMdiff);
% calculate CM (cummulative miosis) and PDRmin (minimum pupil diameter ratio)
    % McLaren, J., Erie, J., & Brubaker, R. (1992). Computerized analysis of pupillograms in studies os alertness. Investigative Ophthalmology & Visual Science, 33(3), 671–676.
CM = sum(1-PDR);
PDRmin = min(PDR);

end

