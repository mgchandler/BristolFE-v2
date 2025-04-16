function s = fn_gaussian_pulse(t, centre_freq, no_cycles)
%SUMMARY
%   Returns a Gaussian pulse with specifed centre frequency and number of
%   cycles (based on -40dB points) given a specified time axis.

sz= size(t);
t = t(:);
beta = 0.01;
ct = mean(t);
T = no_cycles / (2 * centre_freq * sqrt(log(1/beta)));
env = exp(-((t - ct)/T) .^ 2);
i = min(find(env > 0.01));
s = env .* sin(2*pi * centre_freq * (t - ct));
s = [s(i:end); zeros(i - 1, 1)];
s = reshape(s, sz);
end

% This MATLAB code generates a Gaussian pulse with a specified center frequency and number of cycles. Here's a breakdown of each part:
% 
% Input Arguments:
% 
% t: Time axis (1D array).
% centre_freq: Center frequency of the Gaussian pulse.
% no_cycles: Number of cycles of the Gaussian pulse.
% Variable Initialization:
% 
% sz= size(t);: Stores the size of the input time axis.
% t = t(:);: Reshapes the time axis into a column vector for calculation consistency.
% beta = 0.01;: A constant value used in Gaussian pulse calculation.
% ct = mean(t);: Calculates the mean of the time axis, which is used as the center time.
% T = no_cycles / (2 * centre_freq * sqrt(log(1/beta)));: Calculates the width parameter of the Gaussian pulse based on the number of cycles and center frequency.
% Gaussian Envelope Calculation:
% 
% env = exp(-((t - ct)/T) .^ 2);: Computes the Gaussian envelope of the pulse.
% Finding the Starting Point:
% 
% i = min(find(env > 0.01));: Finds the index where the envelope is greater than 0.01, indicating the start of the significant part of the pulse.
% Generating the Pulse:
% 
% s = env .* sin(2*pi * centre_freq * (t - ct));: Generates the Gaussian pulse by modulating the sine wave with the Gaussian envelope.
% Trimming and Reshaping:
% 
% s = [s(i:end); zeros(i - 1, 1)];: Trims the pulse to remove any leading zeros before the significant part of the pulse.
% s = reshape(s, sz);: Reshapes the pulse to match the size of the input time axis.
% Output:
% 
% s: Returns the Gaussian pulse.
% Overall, this code calculates a Gaussian pulse with specified properties and returns it as a function of time.