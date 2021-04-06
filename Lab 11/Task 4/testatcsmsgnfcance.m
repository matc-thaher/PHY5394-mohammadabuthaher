%% Calculate GLRT for Quadratic chirp signal 
% Generalized Likelihood ratio test (GLRT) for a quadratic chirp when only
% the amplitude is unknown.

%% Parameters for data realization
% Number of samples and sampling frequency.
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

%% Supply PSD values
% noise psd using function handle
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);

%% Generate  data realization
% Signal parameters 
a1=10;
a2=3;
a3=3;
A=10;

% Load Data
qcData1 = load("data1.txt");
qcData2 = load("data2.txt");
qcData3 = load("data3.txt");

% data realization for noise and data
noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
dataVec1 = noiseVec+qcData1';
dataVec2 = noiseVec+qcData2';
dataVec3 = noiseVec+qcData3';

%% Plots
% Plot of Data and quadratic chirp
%figure;
%subplot(3,1,1)
%plot(timeVec,dataVec1);
%hold on;
%plot(timeVec,qcData1');
%xlabel('Time (sec)');
%ylabel('Data');
%title('Data realization for calculation of LR');
%subtitle('Data ( Noise + Signal) and First Quadratic Chirp Signal');
%subplot(3,1,2)
%plot(timeVec,dataVec2);
%hold on;
%plot(timeVec,qcData2');
%xlabel('Time (sec)');
%ylabel('Data');
%subtitle('Data ( Noise + Signal) and Second Quadratic Chirp Signal');
%subplot(3,1,3)
%plot(timeVec,dataVec3);
%hold on;
%plot(timeVec,qcData3');
%xlabel('Time (sec)');
%ylabel('Data');
%subtitle('Data ( Noise + Signal) and Third Quadratic Chirp Signal');

% Periodogram of data and signal
%figure;
%kNyq = floor(nSamples/2)+1;
%dataLen = nSamples/sampFreq;
%posFreq = (0:(kNyq-1))*(1/dataLen);
%subplot(3,1,1)
%datFFT1 = abs(fft(dataVec1));
%sigFFT1 = abs(fft(qcData1));
%plot(posFreq,datFFT1(1:kNyq));
%hold on;
%plot(posFreq,sigFFT1(1:kNyq));
%xlabel('Frequency (Hz)');
%ylabel('Periodogram');
%title('Periodogram of data and first quadratic chirp signal');
%subplot(3,1,2)
%datFFT2 = abs(fft(dataVec2));
%sigFFT2 = abs(fft(qcData2));
%plot(posFreq,datFFT2(1:kNyq));
%hold on;
%plot(posFreq,sigFFT2(1:kNyq));
%xlabel('Frequency (Hz)');
%ylabel('Periodogram');
%title('Periodogram of data and second quadratic chirp signal');
%subplot(3,1,3)
%datFFT3 = abs(fft(dataVec3));
%sigFFT3 = abs(fft(qcData3));
%plot(posFreq,datFFT3(1:kNyq));
%hold on;
%plot(posFreq,sigFFT3(1:kNyq));
%xlabel('Frequency (Hz)');
%ylabel('Periodogram');
%title('Periodogram of data and third quadratic chirp signal');

% Spectrogram of data
%figure;
%subplot(3,1,1)
%[S,F,T] = spectrogram(dataVec1,64,60,[],sampFreq);
%imagesc(T,F,abs(S)); axis xy;
%xlabel('Time (sec)');
%ylabel('Frequency (Hz)');
%title('Spectrogram of data (noise + first quadratic chirp signal)');
%subplot(3,1,2)
%[S,F,T] = spectrogram(dataVec2,64,60,[],sampFreq);
%imagesc(T,F,abs(S)); axis xy;
%xlabel('Time (sec)');
%ylabel('Frequency (Hz)');
%title('Spectrogram of data (noise + second quadratic chirp signal)');
%subplot(3,1,3)
%[S,F,T] = spectrogram(dataVec3,64,60,[],sampFreq);
%imagesc(T,F,abs(S)); axis xy;
%xlabel('Time (sec)');
%ylabel('Frequency (Hz)');
%title('Spectrogram of data ( noise + third quadratic chirp signal)');

%% Compute GLRT
% Call the function to compute glrt
glrt1 = glrtqcsig(timeVec, A, sampFreq, dataVec1, psdPosFreq, [a1, a2, a3]);
glrt2 = glrtqcsig(timeVec, A, sampFreq, dataVec2, psdPosFreq, [a1, a2, a3]);
glrt3 = glrtqcsig(timeVec, A, sampFreq, dataVec3, psdPosFreq, [a1, a2, a3]);

fprintf("GLRT value for first QC signal = %d\n", glrt1)
fprintf("GLRT value for second QC signal = %d\n", glrt2)
fprintf("GLRT value for third QC signal = %d\n", glrt3)
%% Likelihood test - null hypothesis
% data realization
noiseRlztn = 2 + 20 .* randn(1, 80000);

% Number of gamma >= gamma(observed)
gamma1 = numel(find(noiseRlztn >= glrt1));
gamma2 = numel(find(noiseRlztn >= glrt2));
gamma3 = numel(find(noiseRlztn >= glrt3));

% probability
pr1 = gamma1/length(noiseRlztn);
pr2 = gamma2/length(noiseRlztn);
pr3 = gamma3/length(noiseRlztn);

fprintf("Estimation significance for first QC signal = %d\n", pr1)
fprintf("Estimation significance for second QC signal = %d\n", pr2)
fprintf("Estimation significance for third QC signal = %d\n", pr3)