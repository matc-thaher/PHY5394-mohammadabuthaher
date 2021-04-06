%% Normalize the linear Transient Chirp signal for a given SNR

% This is the target SNR for the LR
snr = 10;

%% Data generation parameters
nSamples = 2048;
sampFreq = 1024;


%% Generate the signal that is to be normalized
% Signal parameters
A = 5;
t_a = 0.4;
L = 1.3;
f_0 = 2;
f_1 = 5;
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
