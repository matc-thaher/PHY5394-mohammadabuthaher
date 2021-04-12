%% Normalize the linear Transient Chirp signal for a given SNR
%SDM **************************
%Path to signal generation code (should be to lab 2 actually?)
addpath '../Task 1'
% Path added to previous lab containing the iLIGO PSD file
addpath '../../Lab 8/Part 2'
%Replace the path below with that of your own copy of DATASCIENCE_COURSE
addpath /Users/Soumya/Documents/TEMP/DATASCIENCE_COURSE/DETEST/
%******************************
% This is the target SNR for the LR
snr = 10;

%% Data generation parameters
nSamples = 2048;
%SDM**************************
%Since we want to have the LIGO PSD up to at least 700 Hz, the Nyquist rate
%has to be higher. Hence, the sampling frequency has to be higher.
sampFreq = 2048;
%******************************


%% Generate the signal that is to be normalized
% Signal parameters
A = 5;
t_a = 0.2;
L = 1.3;
f_0 = 50;
f_1 = 100;
phase = 0;

% Time samples
timeVec = (0:nSamples-1)/sampFreq;

% Amplitude value does not matter as it will be changed in the normalization
sigVec = atcsmgenltcsig(timeVec,[t_a, t_a + L], A,[f_0,f_1], phase);

%% NOise PSD
% importing the data and modify it
data = load("iLIGOSensitivity.txt", 'ascii');
data = [data; [0, var(data(:,2))]];
data = sortrows(data, 'ascend');
noisePSD = data(:, 2);
freqVal = data(:,1);

%% Generation of the PSD vector for positive DFT values 
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = interp1(freqVal, noisePSD, posFreq);

% Data modification
LIGOdata = table(posFreq', psdPosFreq');
idx = find(LIGOdata.Var1 == 50);
for i = 1:idx
    LIGOdata.Var2(i) = LIGOdata.Var2(idx);
end
PSDVals = LIGOdata.Var2';

%SDM**************************
%Also need to flatten above 700 Hz.
idx = find(posFreq >= 700, 1 );
PSDVals(idx:end)=PSDVals(idx);
figure;
loglog(posFreq,PSDVals);
% PSD sent to innerprodpsd, statgaussnoisegen; not sqrt(PSD) given in
% iLIGOsensitivity.txt
PSDVals = PSDVals.^2;
%*****************************

%% Calculation of the norm
% Norm of signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,sampFreq,PSDVals);
% Normalize signal to specified SNR
sigVec = snr*sigVec/sqrt(normSigSqrd);

%% Test
%Obtain LLR values for multiple noise realizations
nH0Data = 1000;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),PSDVals(:)],100,sampFreq);
    llrH0(lp) = innerprodpsd(noiseVec,sigVec,sampFreq,PSDVals);
end
% figure;
% pwelch(noiseVec,128,[],[],sampFreq);

%Obtain LLR for multiple data (=signal+noise) realizations
nH1Data = 1000;
llrH1 = zeros(1,nH1Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),PSDVals(:)],100,sampFreq);
    % Add normalized signal
    dataVec = noiseVec + sigVec;
    llrH1(lp) = innerprodpsd(dataVec,sigVec,sampFreq,PSDVals);
end
%% SNR estimate
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
