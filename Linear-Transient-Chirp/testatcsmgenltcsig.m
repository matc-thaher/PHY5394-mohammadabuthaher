%% Plot the linear transient chirp signal
% Signal parameters
A = 5;
t_a = 0.7;
L = 3;
f_0 = 2;
f_1 = 5;
phase = 0;
% Maximum instantaneous frequency after t_a + L sec is
maxinstFreq = f_0 + 2 * f_1 * L;
samplFreq = 5 * maxinstFreq;
samplIntrvl = 1/samplFreq;

% Time samples
timedata = -5:samplIntrvl:5;
% Number of samples
nSamples = length(timedata);
% Generate the signal
sigVec = atcsmgenltcsig(timedata,[t_a, t_a + L], A,[f_0,f_1], phase);
%Plot the signal 
figure;
plot(timedata,sigVec, '* -')
xlabel("Time in sec")
ylabel("s(t)")
title("Linear Transient Chirp Signal")

%Plot the Periodogram
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
