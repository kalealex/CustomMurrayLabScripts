function [ ageInMonths ] = AK_calculateInterviewAge_forNDAR( DateOfBirth, InterviewDate )
%AK_calculateInterviewAgeForNDAR calculates the age (in months) of a
%participant for the purposes of reporting data to NDAR.
%   INPUT:
%       DateOfBirth: date of birth as a formatted date string or a datenum
%       InterviewDate: date of data collection as a formatted date string
%           or a datenum
%   OUTPUT:
%       ageInMonths: age in months at the time of data collection as a
%           double, rounded to the nearest integer according to NDAR
%           standard (round up if days >= 15)

% check input
if nargin < 2
    error('AK_calculateInterviewAgeForNDAR requires two arguments, date of birth and interview date, as either date strings or datenums.');
end

% convert to datenum units
if ischar(DateOfBirth)
    DateOfBirth = datenum(DateOfBirth);
end
if ischar(InterviewDate)
    InterviewDate = datenum(InterviewDate);
end

% calculate age as differences between dates
difference = InterviewDate - DateOfBirth;

% break down difference into years, months, and days
years = str2double(datestr(difference,'yyyy'));
months = str2double(datestr(difference,'mm'));
days = str2double(datestr(difference,'dd'));

% add together, converting years to months and rounding days according to
% NDAR stardard (round up if days >= 15)
ageInMonths = (12 * years) + months + floor(days / 15);

end

