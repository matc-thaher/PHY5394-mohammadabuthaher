%% Normalization of linear transient chirp signal for a given SNR
% changed the path according to the local repository
addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\DATASCIENCE_COURSE\DETEST'
addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\DATASCIENCE_COURSE\NOISE'
% This is the target SNR for the LR
snr = 10;

%%
% Data generation parameters
nSamples = 2048;
sampFreq = 1024;


%% Generate the signal that is to be normalized
% Signal parameters
A = 5;
t_a = 0.4;
L = 1.3;
f_0 = 50;
f_1 = 70;
phase = 0;

% Time samples
timeVec = (0:nSamples-1)/sampFreq;

% Amplitude value does not matter as it will be changed in the normalization
sigVec = atcsmgenltcsig(timeVec,[t_a, t_a + L], A,[f_0,f_1], phase);

%% Noise
% Noise PSD
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;

%% Generation of PSD vector for positive DFT values
% PSD Vector
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);
figure;
plot(posFreq,psdPosFreq);
axis([0,posFreq(end),0,max(psdPosFreq)]);
xlabel('Frequency (Hz)');
ylabel('PSD ((data unit)^2/Hz)');

%% Calculation of the norm
% Norm of signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,sampFreq,psdPosFreq);
% Normalize signal to specified SNR
sigVec = snr*sigVec/sqrt(normSigSqrd);

%% Test
%Obtain LLR values for multiple noise realizations
nH0Data = 1000;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    llrH0(lp) = innerprodpsd(noiseVec,sigVec,sampFreq,psdPosFreq);
end
%Obtain LLR for multiple data (=signal+noise) realizations
nH1Data = 1000;
llrH1 = zeros(1,nH1Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    % Add normalized signal
    dataVec = noiseVec + sigVec;
    llrH1(lp) = innerprodpsd(dataVec,sigVec,sampFreq,psdPosFreq);
end
%%
% Signal to noise ratio estimate
estSNR = (mean(llrH1)-mean(llrH0))/std(llrH0);

figure;
histogram(llrH0);
hold on;
histogram(llrH1);
xlabel('LLR');
ylabel('Counts');
legend('H_0','H_1');
title(['Estimated SNR = ',num2str(estSNR)]);

%% Plot of Data realization and signal
% A data realization
figure;
plot(timeVec,dataVec);
hold on;
plot(timeVec,sigVec);
xlabel('Time (sec)');
ylabel('Data');
title('Data realization and signal')

%% Plot periodogram
% FFT of signal
fftSig = fft(sigVec);
% fft of data
fftDat = fft(dataVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);
fftDat = fftDat(1:kNyq);

% plot
figure;
plot(posFreq, abs(fftDat), 'b -');
hold on
plot(posFreq,abs(fftSig), 'r -');
hold off
xlabel("Frequency(Hz)")
ylabel("Magnitude")
title("Periodogram of Data and Signal")

%% Plot a spectrogram

winLen = 0.034;%sec
ovrlp = 0.030;%sec
%Convert to integer number of samples 
winLenSmpls = floor(2*winLen*sampFreq);
ovrlpSmpls = floor(2*ovrlp*sampFreq);
[S,F,T]=spectrogram(dataVec,winLenSmpls,ovrlpSmpls,[],sampFreq);

% Plot spectrogram
figure;
imagesc(T,F,abs(S)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title("Spectrogram of Data")


