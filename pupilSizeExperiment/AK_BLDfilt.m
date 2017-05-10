function [ filtTimeSeries, velo, accel, jerk ] = AK_BLDfilt( timeSeries, passFreq, stopFreq, sampleRate )
%AK_BLDfilt applies a band-limited differentiating filter to an array of
%time series data, generating low-pass filtered versions of the time series
%input, the velocity of the time series, the acceleration of the time
%series, and the jerk of the time series.
%   INPUT:
%       timeSeries: an array of time series data
%       passFreq: a frequency below which the time series data will be
%           unchanged in the Fourier domain
%       stopFreq: a frequency above which the time series data will be
%           mostly attenuated in the Fourier domain (the filter is applied
%           as a trapezoidal window in Fourier space, with no filtering
%           below the pass band frequency, linear interpolation between the
%           pass band and the stop band frequencies, and ripling above the
%           stop band to mitigate ringing artifacts)
%       sampleRate: the sample rate of the time series data (hz)
%   OUTPUT:
%       filtTimeSeries: a low-pass filtered version of the time series data
%       velo: the velocity (fisrt derivative) of the time series
%       accel: the acceleration (second derivative) of the time series
%       jerk: the jerk (third derivative) of the time series

% 170510: AMK changed output from [velo, filtTimeSeries, jerk] to [filtTimeSeries, velo, accel, jerk]

% check inputs
if nargin < 4
    error(['AK_BLDfilt requires input:'...
        ' 1) a cell array of trials, each containing a trace of X coordinates for eye tracking during ssMotion, or a single array of contiguous time series data;'...
        ' 2) passband frequency as a double;'...
        ' 3) stopband frequency as a double;'...
        ' 4) sampling rate as a double']) 
end
if ~iscell(timeSeries)
    % output should be matrices of time series data
    matOut = true;
    timeSeries = {timeSeries};
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
for iTr = 1:length(timeSeries)
    % differentiation
    velo{iTr} = filter(df,[timeSeries{iTr}; zeros(D,1)]); % trace' filtered
    velo{iTr} = velo{iTr}(D+1:end); % correct delay
    accel{iTr} = filter(df,[velo{iTr}; zeros(D,1)]); % trace'' filtered
    accel{iTr} = accel{iTr}(D+1:end); % correct delay
    velo{iTr} = velo{iTr}/dt; % divide by time differential
    accel{iTr} = accel{iTr}/dt^2;
    % integration
    filtTimeSeries{iTr} = filter(1,[1 -0.999],velo{iTr}); 
%     velocity{iTr} = filter(1,[1 -0.999],accel{iTr});
    filtTimeSeries{iTr} = filtTimeSeries{iTr}*dt; % multiply by time differential
%     velocity{iTr} = velocity{iTr}*dt^2;
    % determine jerk by differentiating unfiltered acceleration
    jerk{iTr} = diff(timeSeries{iTr},3); % trace''' unfiltered
    jerk{iTr} = jerk{iTr}/dt^3; % divide by time differential
end

% keep input and output as the same data type
if matOut
    filtTimeSeries = cell2mat(filtTimeSeries);
    velo = cell2mat(velo);
    accel = cell2mat(accel);
    jerk = cell2mat(jerk);
end

end

