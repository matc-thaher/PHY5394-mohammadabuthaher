%% Whitenning the data
%Sampling frequency for noise realization
sampFreq = 1024; %Hz
%Number of samples to generate
nSamples = 16384;

% Filter order
fltrOrdr = 500;

%% Estimating the data
% Loading the data as table
noiseData = readtable("testData.txt");

% Changing variable names
noiseData.Properties.VariableNames = ["Time", "Strain"];

% time data
timeVec = noiseData.Time;

% Selecting time before 5.0s
mnsData = noiseData(noiseData.Time <= 5, :);

% Welch estimated psd
[psdVecD, freqVecD] = pwelch(mnsData.Strain);


%% Whitenning the noise realization
inNoise = noiseData.Strain;
outNoise = whitenning(inNoise,[freqVecD(:),psdVecD(:)],fltrOrdr,sampFreq);


%% Estimate and plotting the PSD
[psd_in,freq_in]=pwelch(inNoise, 256,[],[],sampFreq);
[psd_out,freq_out]=pwelch(outNoise, 256,[],[],sampFreq);

figure;
subplot(2,1,1)
plot(freq_in, psd_in);
xlabel('Frequency (Hz)');
ylabel('PSD');
title('Welch PSD estimate for Input Noise')

subplot(2,1,2)
plot(freq_out,psd_out);
xlabel('Frequency (Hz)');
ylabel('PSD');
title('Welch PSD estimate for Output Noise')


%% Plot a spectrogram
% Parameters for spectrogram
winLen = 0.25;%sec
ovrlp = 0.2;%sec

% Convert to integer number of samples 
winLenSmpls = floor(2*winLen*sampFreq);
ovrlpSmpls = floor(2*ovrlp*sampFreq);
[Si,Fi,Ti]=spectrogram(inNoise,winLenSmpls,ovrlpSmpls,[],sampFreq);
[So,Fo,To]=spectrogram(outNoise,winLenSmpls,ovrlpSmpls,[],sampFreq);
figure;
subplot(2,1,1)
imagesc(Ti,Fi,abs(Si)); axis xy;
xlabel('Time(s)')
ylabel('Frequency (hz)')
title('Spectrogram of input noise')

subplot(2,1,2)
imagesc(To, Fo, abs(So)); axis xy;
xlabel('Time(s)')
ylabel('Frequency (hz)')
title('Spectrogram of output noise')

%% Plot the time series of noise realization
figure;
subplot(2,1,1)
plot(timeVec,inNoise);
xlabel('Time(s)')
ylabel('Strain')
title('Time series of input noise')

subplot(2,1,2)
plot(timeVec,outNoise);
xlabel('Time(s)')
ylabel('Strain')
title('Time series of output noise')

%% Function for whitenning 
function outNoise = whitenning(inNoise,psdVals,fltrOrdr,sampFreq)
%Whitenning the noise data
%Y = STATGAUSSNOISEGEN(N,PSD,O,Fs)
%Give a whitened data Y of an input noise. Fs is the sampling frequency
%of Y. PSD is a vector containing frequencies and the corresponding
%PSD values in the first and second columns for first 5.0s respectively.
%The frequencies must start from 0 and end at 1. The order of the FIR 
%filter to be used is given by O.

%Mohammad Abu Thaher Chowdhury, Mar 2021

% Design FIR filter with T(f)= square root of PSD
freqVec = psdVals(:,1);
sqrtPSD = sqrt(psdVals(:,2));
rtiosampfreqVec = sampFreq/freqVec(end);
b = fir2(fltrOrdr,((rtiosampfreqVec * freqVec)/sampFreq),sqrtPSD);

% Passing noise data through the designed filter
rng('default'); 
outNoise = sqrt(sampFreq)*fftfilt(b,inNoise);
end