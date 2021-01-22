%% Plot the quadratic chirp signal
% Signal parameters
A = 10;
t_a = 0.7;
L = 10;
f_0 = 5;
f_1 = 7;
phase = 0;
% Instantaneous frequency after 1 sec is
instFreq = f_0 + 2 * f_1 * (1 - t_a);
samplFreq = instFreq;
samplIntrvl = 1/samplFreq;

% Time samples
timedata = t_a - 2.7:samplIntrvl: t_a + L + 5;
% Number of samples
nSamples = length(timedata);
% Generate the signal
sigVec = test_4(timedata,[t_a, t_a + L], A,[f_0,f_1], phase);
%Plot the signal 
figure;
plot(timedata,sigVec, '* -')
%Plot the periodogram
%--------------
%Length of data 
dataLen = timedata(end)-timedata(1);
%DFT sample corresponding to Nyquist frequency
kNyq = floor(nSamples/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/dataLen);
% FFT of signal
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

%Plot periodogram
figure;
plot(posFreq,abs(fftSig));
