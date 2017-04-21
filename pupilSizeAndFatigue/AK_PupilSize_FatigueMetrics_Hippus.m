function [ PowerLowFreq ] = AK_PupilSize_FatigueMetrics_Hippus( pupilSize, sampleRate )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% determine length of signal; should be even
L = length(pupilSize);
% pad fft to nearest power of 2
p = nextpow2(L);
n = 2^p;


% compute the Fourier transform of the signal.
Y = fft(pupilSize,n);
% compute the two-sided spectrum P2; then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
% define the frequency domain f 
f = sampleRate*(0:(L/2))/L;
% integrate with respect to f within the limits of integration 0 - 0.8 hz
        % Lüdtke, H., Wilhelm, B., Adler, M., Schaeffel, F., & Wilhelm, H. (1998). Mathematical procedures in data recording and processing of pupillary fatigue waves. Vision Research, 38(19), 2889–2896. http://doi.org/10.1016/S0042-6989(98)00081-9
[~,limLidx] = min(abs(f-0)); % idx for limits of integration
[~,limUidx] = min(abs(f-0.8));
PowerLowFreq = trapz(P1(limLidx:limUidx));

end

