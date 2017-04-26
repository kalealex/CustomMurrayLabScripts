function AK_prepareNDARsubmission( submission_dir, NDAR_subinfo_xlsfile )
%AK_prepareNDARsubmission is a function which runs a set of other functions
%which prepare different parts of the NDAR submission.
%   INPUT:
%       submission_dir: (string) should be the full path name of a folder where the submission materials are to be saved
%       NDAR_subinfo_xlsfile (string) should be the full path name of a .xls file containing the columns: SUB ID, GUID, Date of Birth, Gender
%           this information can be found on REDCap under the report
%           labeled Subinfo_NDAR for the GABA project (this information is
%           typically added to REDCap by the Bernier Lab coordinators)
%   FOLLOW UP:
%       After running this function, you will need to make a few manual
%       changes to the .csv files it creates in the designated directory.
%       Open each file and use the find/replace tool to replace the symbols
%        ' and [] with nothing. This will simply remove the Matlab
%        formatting for strings and empty arrays. Next, for files with bad 
%        .prt file links, the fields 'experiment_id', 'image_extent4' and 
%        'study' will need to be filled in manually. Since these are all 
%        older subjects (i.e., G102, G105, G307, G311), 
%        the easiest way to do this is to cut and paste these columns from 
%        past submissions at the rows that correspond to these specific 
%        subjects. The validation tool will catch these errors if these 
%        empty cells are not properly filled. Lastly, you will need to use
%        the validation tool, available for download online at 
%       https://ndar.nih.gov/tools_validation_tool.html, to check the 
%       validity of the submission folder and then upload it to NDAR.
%       This will complete the submission process.


% check inputs
if nargin < 2 || nargin > 2
    error(['AK_prepareNDARsubmission requires two input arguments:',...
        'submission_dir (string) should be the full path name of a folder where the submission materials are to be saved;'...
        'NDAR_subinfo_xlsfile (string) should be the full path name of a .xls file containing the columns: SUB ID, GUID, Age (months), Gender']);
end

% copy the appropriate files to the submission folder
AK_copyMRIdata_forNDAR(submission_dir, NDAR_subinfo_xlsfile)

% prepare .xls file containing fMRI meta-data
AK_fMRImetadata_forNDAR(fullfile(submission_dir,'image03.csv'), NDAR_subinfo_xlsfile)

% prepare .xls file containing MRS meta-data
AK_MRSmetadata_forNDAR(fullfile(submission_dir,'cu_spectroscopy01.csv'), NDAR_subinfo_xlsfile)

% prepare .xls file containing psychophysics data
AK_PsychophysicsData_forNDAR(fullfile(submission_dir,'psef01.csv'), NDAR_subinfo_xlsfile)

end

