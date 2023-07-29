clear all;

t = linspace(0,1,100); % time vector
% Read and populate noise and tap coefficient vectors
coeffs = csvread('NewFilters.csv');
Noise2 = csvread('noiseVector2.csv');
Noise3 = csvread('noiseVector3.csv');
Noise4 = csvread('noiseVector4.csv');
OriginalNoise = csvread('noiseVector.csv');

% Create discrete ideal input signal vector
cleanSignal = sin(2*pi*2*t);
%plot(t,cleanSignal)
% Create noisy version of input signal vector
noisySignal = cleanSignal + 0.1 * OriginalNoise;
noisySignal2 = cleanSignal + 0.1 * Noise2;
noisySignal3 = cleanSignal + 0.1 * Noise3;
noisySignal4 = cleanSignal + 0.1 * Noise4;
%plot(t,noisySignal)
% Create and fill matrix with ASE for each filter/signal pair
ASE_matrix = zeros(5,3)
for k = 1:5
 ASE_matrix(k,1) = ComputeASE(cleanSignal, noisySignal2, coeffs(k,:));
 ASE_matrix(k,2) = ComputeASE(cleanSignal, noisySignal3, coeffs(k,:));
 ASE_matrix(k,3) = ComputeASE(cleanSignal, noisySignal4, coeffs(k,:));
end
% display the matrix
ASE_matrix


% Plot original and filtered signals for each noise set
% using set 1 of filter coefficients
[A, OS] = ComputeASE(cleanSignal, noisySignal, coeffs(1,:))
subplot(2,4,1)
plot(t(5:end), noisySignal(5:end))
title('Noisy Training Signal')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on
subplot(2,4,5)
plot(t(5:end), OS)
title('Filtered Training Signal')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on
[A, OS] = ComputeASE(cleanSignal, noisySignal2, coeffs(1,:))
subplot(2,4,2)
plot(t(5:end), noisySignal2(5:end))
title('Noisy Test Signal 1')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on
subplot(2,4,6)
plot(t(5:end), OS)
title('Filtered Test Signal 1')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on
[A, OS] = ComputeASE(cleanSignal, noisySignal3, coeffs(1,:))
subplot(2,4,3)
plot(t(5:end), noisySignal3(5:end))
title('Noisy Test Signal 2')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on
subplot(2,4,7)
plot(t(5:end), OS)
title('Filtered Test Signal 2')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on
[A, OS] = ComputeASE(cleanSignal, noisySignal4, coeffs(1,:))
subplot(2,4,4)
plot(t(5:end), noisySignal4(5:end))
title('Noisy Test Signal 3')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on
subplot(2,4,8)
plot(t(5:end), OS)
title('Filtered Test Signal 3')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on
% This function passes the noisy signal through the desired filter,
% compares the output to the original clean signal, and computes the ASE
% This function is the same as the function from the original
% program except for the required transposition of the tap
% coefficient matrix
function [ASE, outputSignal] = ComputeASE(cleanSignal, noisySignal,
coefficients)

 iterations = length(noisySignal) - length(coefficients) + 1;
 outputSignal = zeros(1, iterations);

 for i=1:iterations
 outputSignal(i) = flip(noisySignal(i:(i + length(coefficients) - 1)))
...
 * transpose(coefficients);
 end

 errorVector = cleanSignal(length(coefficients):end) - outputSignal;
 ASE = sum(errorVector .^ 2) / length(errorVector);

end