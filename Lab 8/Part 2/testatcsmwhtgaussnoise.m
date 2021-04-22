%% Generate a WGN realization and filter it using
%Sampling frequency for noise realization
sampFreq = 4096; %Hz
%SDM: Number of samples
nSamples = 20*sampFreq;

% importing the data
ligoData = readtable('iLIGOSensitivity.txt');
ligoData.Properties.VariableNames = ["Frequency", "sqrtPSD"];

%% Data Modifications
% Adding frequency 0, 50, fs/2 and 700, along with psd Value, to the data by linear interpolation
freqVal = 0:(sampFreq/2);
PSDval = interp1(ligoData.Frequency, ligoData.sqrtPSD, freqVal);
T2 = table(freqVal, PSDval);
idx3 = find(T2.freqVal == 0 | T2.freqVal == sampFreq/2| T2.freqVal == 50 | T2.freqVal == 700);
Frequency = (T2.freqVal(idx3))';
sqrtPSD = (T2.PSDval(idx3))';
T3 = table(Frequency, sqrtPSD);
ligoData = [ligoData;T3];
ligoData = sortrows(ligoData,'Frequency','ascend');

% Making the values constant before and after 50 hz and 700 hz respectively
for idx1 = find(ligoData.Frequency <= 50)
ligoData.sqrtPSD(idx1) = ligoData.sqrtPSD(idx1(end));
end
for idx2 = find(ligoData.Frequency >= 700)
ligoData.sqrtPSD(idx2) = ligoData.sqrtPSD(idx2(1));
end


% checking constant data before 50 hz and 700 hz with logarithmic plot
figure;
loglog(ligoData.Frequency,ligoData.sqrtPSD)
xlabel("Frequency ")
ylabel("Square root of (s_n(f)")
title("Data Modification")
ylim([(1 * 10^-23) (16 * 10^-22)])

%% SDM: First limit the frequency values
SDM_fVals = ligoData.Frequency;
SDM_fVals(SDM_fVals > sampFreq/2) = [];
SDM_sqrtPSD = ligoData.sqrtPSD;
SDM_sqrtPSD = SDM_sqrtPSD(1:length(SDM_fVals));

%% Designing a filter using PSD data
fltrOrdr = 300;
% ratioVal = sampFreq/ligoData.Frequency(end);
% b = fir2(fltrOrdr,ligoData.Frequency * ratioVal/(sampFreq),ligoData.sqrtPSD);
b = fir2(fltrOrdr,SDM_fVals/(sampFreq/2),SDM_sqrtPSD);
%% Generating White Gaussian Noise from PSD values
% Creating noise by calling the function
% noise1 = genNoisefrmPSD(ligoData.sqrtPSD);

% concatenating noise tables into one noise table
% inNoise = vertcat(noise1, noise2, noise3, noise4, noise5);
inNoise = randn(1,nSamples);

% Passed the noise through the filter, built by using PSD value
outNoise = fftfilt(b,inNoise);

figure;
plot((0:(nSamples-1))/sampFreq, outNoise);

%% Estimate the PSD
[pxxo, fo]= pwelch(outNoise, sampFreq, [], [], sampFreq);
figure;
%pwelch(outNoise)
loglog(fo,pxxo);

%% Function for creating white gaussian noise from PSD
function testnoise = genNoisefrmPSD(a)
% N = GENNOISEFRMPSD(A)
% Generate a noise realization using PSD data. Here, A represents
% square root of power spectral density.

%Mohammad Abu Thaher Chowdhury, March 2021

sData = randn(length(a), 1);  %generating noise code (outside function) examples in comment
sDft = abs(fft(sData));       %sData4 = randn(length(ligoData.sqrtPSD), 1);
sMult = sDft .* a;            %s4dft = abs(fft(sData4));
testnoise = abs(ifft(sMult)); %s4mult = s4dft * ligoData.sqrtPSD;
                              %signal4 = abs(ifft(s4mult));
end

