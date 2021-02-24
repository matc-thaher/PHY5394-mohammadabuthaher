%% Plot the spectrogram, periodogram and the linear transient chirp signal
% Signal parameters
A = 2;
t_a = 0.4;
L = 1.2;
f_0 = 5;
f_1 = 15;
phase = 0;
% Instantaneous frequency after 1 sec is
instFreq = f_0 + 2 * f_1 * L;
samplFreq = 5 * instFreq;
samplIntrvl = 1/samplFreq;

% Time samples
timedata = 0:samplIntrvl: 2;
% Number of samples
nSamples = length(timedata);
% Generate the signal
sigVec = atcsmsptgmltcsig(timedata,[t_a, t_a + L], A,[f_0,f_1], phase);
%Plot the signal 
figure;
plot(timedata,sigVec, '* -')
title("Linear Transient Chirp Signal")
xlabel("Time(s)")
ylabel("s(t)")

%% Plot the periodogram
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
plot(posFreq,abs(fftSig), 'm -');
xlabel("Frequency(Hz) (only positive values)")
ylabel("Magnitude of s(t)")
title("Fourier Transform of Linear Transient Chirp Signal")

%% Plot a spectrogram
%----------------
winLen = 0.25;%sec
ovrlp = 0.2;%sec
%Convert to integer number of samples 
winLenSmpls = floor(2*winLen*samplFreq);
ovrlpSmpls = floor(2*ovrlp*samplFreq);
[S,F,T]=spectrogram(sigVec,winLenSmpls,ovrlpSmpls,[],samplFreq);
figure;
imagesc(T,F,abs(S)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title("Spectrogram of Linear Transient Chirp Signal")
