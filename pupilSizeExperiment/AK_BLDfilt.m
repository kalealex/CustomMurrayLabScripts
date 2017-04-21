function [ velo, filtTrace, jerk ] = AK_BLDfilt( X_trace, passFreq, stopFreq, sampleRate )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% check inputs
if nargin < 4
    error(['AK_BLDfilt requires input:'...
        ' 1) a cell array of trials, each containing a trace of X coordinates for eye tracking during ssMotion, or a single array of contiguous time series data;'...
        ' 2) passband frequency as a double;'...
        ' 3) stopband frequency as a double;'...
        ' 4) sampling rate as a double']) 
end
if ~iscell(X_trace)
    % output should be matrices of time series data
    matOut = true;
    X_trace = {X_trace};
else
    matOut = false;
end

% setup filter
Fs  = sampleRate; % sample rate
dt = 1/Fs; % time differential
df = designfilt('differentiatorfir', ... % Response type
       'FilterOrder',50, ...            % Filter order
       'PassbandFrequency',passFreq, ...     % Frequency constraints
       'StopbandFrequency',stopFreq, ...
       'DesignMethod','equiripple', ... % Design method
       'SampleRate',Fs);               % Sample rate
D = mean(grpdelay(df)); % filter delay

% apply differentiating filter
for iTr = 1:length(X_trace)
    % differentiation
    velo{iTr} = filter(df,[X_trace{iTr}; zeros(D,1)]); % trace' filtered
    velo{iTr} = velo{iTr}(D+1:end); % correct delay
    accel{iTr} = filter(df,[velo{iTr}; zeros(D,1)]); % trace'' filtered
    accel{iTr} = accel{iTr}(D+1:end); % correct delay
    velo{iTr} = velo{iTr}/dt; % divide by time differential
    accel{iTr} = accel{iTr}/dt^2;
    % integration
    filtTrace{iTr} = filter(1,[1 -0.999],velo{iTr}); 
    velocity{iTr} = filter(1,[1 -0.999],accel{iTr});
    filtTrace{iTr} = filtTrace{iTr}*dt; % multiply by time differential
    velocity{iTr} = velocity{iTr}*dt^2;
    % determine jerk by differentiating unfiltered acceleration
    jerk{iTr} = diff(X_trace{iTr},3); % trace''' unfiltered
    jerk{iTr} = jerk{iTr}/dt^3; % divide by time differential
end

% keep input and output as the same data type
if matOut
    velo = cell2mat(velo);
    filtTrace = cell2mat(filtTrace);
    jerk = cell2mat(jerk);
end

end

