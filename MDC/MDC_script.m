% Adding path for the specific functions and data
 addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\DATASCIENCE_COURSE\MDC'
 addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\PHY5394-mohammadabuthaher\Final labs\Lab 4'
 addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\DATASCIENCE_COURSE\DETEST'
 addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\DATASCIENCE_COURSE\SIGNALS'
 addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\DATASCIENCE_COURSE\NOISE'
 addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\SDMBIGDAT19\CODES'
 addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\PHY5394-mohammadabuthaher\Lab 11\Task 4'

% Loading the data
trData = load('TrainingData.mat');
anData = load('analysisData.mat');

%% Parameters
nSamples = length(anData.dataVec);
% Search range of phase coefficients
rmin = [40, 1, 1];
rmax = [100, 50, 15];

% Number of independent PSO runs
nRuns = 8;


%% Whitenning the noise realization
%SDM***********************
%You estimated the whitening filter from training data but did not whiten
%the noise but estimated the whitened PSD from the data containing the
%signal. This is messed up: The final PSD in the inner product is the
%estimated PSD of *noise*.
%I have removed the whitening part below.
% filter
% fltrOrdr = 100;
% [psdVecD, freqVecD] = pwelch(trData.trainData, nSamples, [], [], nSamples);
% sqrtPSD = sqrt(psdVecD);
% ratio = (anData.sampFreq/2)/freqVecD(end);
% b = fir2(fltrOrdr,((ratio * freqVecD)/(anData.sampFreq/2)),sqrtPSD);
% 
% whitening the data
inData = anData.dataVec;
%outData = fftfilt(b,inData);
outData = inData;

% estimating psd
%The estimated PSD from noise-only data because it is supposed to be the
%PSD of *noise*. Because you estimated the PSD from the shorter analysis
%data, you got a very noisy estimate that threw the parameter estimate off.
%[psdOut, freqOut] = pwelch(outData,nSamples,[],[], anData.sampFreq);
[psdOut, freqOut] = pwelch(trData.trainData,nSamples,[],[], anData.sampFreq);
% dataLen = nSamples/anData.sampFreq;
% kNyq = floor(nSamples/2)+1;
% posFreq = (0:(kNyq-1))*(1/dataLen);
posFreq = freqOut;
psdPosFreq = psdOut';
%***************************

%% 
% time vector
timeVec = (0:(nSamples-1))/nSamples;

% Input parameters for GLRTQCSIGPSO
inParams = struct('dataX', timeVec,...
                  'dataY', outData,...
                  'dataXSq',timeVec.^2,...
                  'dataXCb',timeVec.^3,...
                  'psdVal', psdPosFreq,...
                  'sampFreq',anData.sampFreq,...
                  'rmin',rmin,...
                  'rmax',rmax);
              
% GLRTQCSIGPSO runs PSO on the GLRTQCSIG4PSO4 fitness function. As an
% illustration of usage, we change one of the PSO parameters from its
% default value.
outStruct = glrtqcsigpso(inParams,struct('maxSteps',2000),nRuns);

%% Computing glrt
glrtPSO = -outStruct.bestFitness;

%mle
llr = sqrt(glrtPSO);

%% Estimation
% data realization parameters
m = 2000;
glrtH0PSO = zeros(1,m);

for s = 1:m
outNoise = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,anData.sampFreq);
%SDM****************
%This is incorrect: The noise realizations should be processed through
%exactly the same method as was used on the analysis data vector. This
%means running PSO on each noise realization. This means that you don't need to use
%the estimated signal parameters.
glrt4work= glrtqcsig(timeVec, 10, anData.sampFreq, outNoise, psdPosFreq, outStruct.bestQcCoefs);
%*******************
glrtH0PSO(s) = glrt4work;
end


% Number of gamma >= gamma(observed)
gamma = numel(find(glrtH0PSO >= glrtPSO));

% estimating significance, signal exists as pr < 0.5
pr = gamma/length(glrtH0PSO);

 %% plot
 
figure;
hold on;
plot(timeVec,inParams.dataY,'.');
plot(timeVec,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
xlabel('Time in second')
title('Estimation of Signal')
legend('Data',...
       ['Estimated signal: ',num2str(nRuns),' runs']);

disp(['Estimated parameters: a1=',num2str(outStruct.bestQcCoefs(1)),...
                             '; a2=',num2str(outStruct.bestQcCoefs(2)),...
                             '; a3=',num2str(outStruct.bestQcCoefs(3)),...
                             '; Estimated significance: ',num2str(pr)])
                         
% Spectrogram of data before filer and after filter
% figure;
% subplot(2,1,1)
% [S,F,T] = spectrogram(inData,64,60,[],anData.sampFreq);
% imagesc(T,F,abs(S)); axis xy;
% xlabel('Time (sec)');
% ylabel('Frequency (Hz)');
% title('Spectrogram of input data');
% subplot(2,1,2)
% [S,F,T] = spectrogram(outData,64,60,[],anData.sampFreq);
% imagesc(T,F,abs(S)); axis xy;
% xlabel('Time (sec)');
% ylabel('Frequency (Hz)');
% title('Spectrogram of output data');