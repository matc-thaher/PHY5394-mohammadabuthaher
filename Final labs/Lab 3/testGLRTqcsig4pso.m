%% Calculate GLRT for Quadratic chirp signal 
% Generalized Likelihood ratio test (GLRT) for a quadratic chirp when only
% the amplitude is unknown.
addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\DATASCIENCE_COURSE\DETEST'
addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\DATASCIENCE_COURSE\NOISE'
addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\DATASCIENCE_COURSE\SIGNALS'
addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\SDMBIGDAT19\CODES'
addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\PHY5394-mohammadabuthaher\Lab 11\Task 3'
%% Parameters for data realization
% Number of samples and sampling frequency.
nSamples = 2048;
sampFreq = 2048;
timeVec = (0:(nSamples-1))/sampFreq;

%% Supply PSD values
% This is the noise psd we will use.
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
snr=10;

% Generating QC Data
qc = crcbgenqcsig(timeVec, snr, [a1,a2,a3]);
[qcData,~] = normsig4psd(qc,sampFreq,psdPosFreq,snr);

% data realization
noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
dataVec1 = noiseVec+qcData;

% Array for xVec
A = linspace(1,15, length(timeVec));
A2 = [1 5];
A3 = [1 5];
x1 = ones(1, length(A));
x2 = ones(1, length(A));
x3 = ones(1, length(A));
x1(1,:) = (A - min(A))./ (max(A) - min(A));
x2(1,:) = (a2 - min(A2))./ (max(A2) - min(A2));
x3(1,:) = (a3 - min(A3))./ (max(A3) - min(A3));
xVec = [x1;x2;x3]';

% struct of parameters
ffparams = struct('dataY', dataVec1,...
    'dataX', timeVec,...
    'dataXSq', timeVec.^2,...
    'dataXCb', timeVec.^3,...
    'psdVal', psdPosFreq,...
    'snr', snr,...
    'sampFreq', sampFreq,...
    'rmin', [1,1,1],...
    'rmax', [15,5,5]); 
%% Fitness Value
% Compute the fitness value
fitVal = glrtqcsig4pso(xVec, ffparams);

% Plot of Fitness value
figure;
plot(A,fitVal, '. -')
xlabel('A_i (Array of x^i_1)')
ylabel('Fitness Value')
title('Fitness value vs A_i')