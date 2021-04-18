%% Plot the linear transient chirp signal
% Signal parameters
A = 10;
%t_a = 0.2;
%L = 0.6;
%f_0 = 50;
%f_1 = 70;
%phase = 0;

% parameters
nSamples = 2048;
samplFreq = 2048;
% Time samples
timedata = (0:nSamples-1)/samplFreq;

% struct of the vectors
%sigParameters = struct('sTsig', t_a, 'fTsig', t_a + L, 'inFreq', f_0, 'fnFreq', f_1, 'phase', 0);
sigParameters = struct('sTsig', 0.2, 'fTsig', 0.8, 'inFreq', 50, 'fnFreq', 70, 'phase', 0);
% Generate the signal
sigVec = atcsmgenltcsignew(timedata,A, sigParameters);

%% Plots
% Plot the time series 
figure;
plot(timedata,sigVec, '* -')
xlabel("Time in sec")
ylabel("s(t)")
title("Linear Transient Chirp Signal")

% Parameters for the Periodogram
% Length of data 
dataLen = timedata(end)-timedata(1);
% DFT sample corresponding to Nyquist frequency
kNyq = floor(nSamples/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/dataLen);
% FFT of signal
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

%Plot periodogram
figure;
plot(posFreq,abs(fftSig), 'm -');
xlabel("Frequency in Hz(only positive value)")
ylabel("Magnitude of signal")
title("Periodogram of Linear Transient Chirp Signal")
