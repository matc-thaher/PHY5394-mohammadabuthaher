%% Minimize the fitness function CRCBQCFITFUNC using PSO
%addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\PHY5394-mohammadabuthaher\Lab 11\Task 3'
%addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\SDMBIGDAT19\CODES'
%addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\PHY5394-mohammadabuthaher\Final labs\Lab 3'
% Data length
nSamples = 512;
% Sampling frequency
sampFreq = 512;
% Signal to noise ratio of the true signal
snr = 10;
% Phase coefficients parameters of the true signal
a1 = 10;
a2 = 3;
a3 = 3;

% Search range of phase coefficients
rmin = [1, 1, 1];
rmax = [180, 10, 10];

% Number of independent PSO runs
nRuns = 8;

%% Supply PSD values
% This is the noise psd we will use.
noisePSD = @(f) (f>=50 & f<=100).*(f-50).*(100-f)/625 + 1;
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);
%% Do not change below
% Generate data realization
timeVec = (0:(nSamples-1))/sampFreq;
% Reset random number generator to generate the same noise realization,
% otherwise comment this line out
rng('default');

% Generate data realization

% Generating QC Data
qc = crcbgenqcsig(timeVec, snr, [a1,a2,a3]);
[qcData,~] = normsig4psd(qc,sampFreq,psdPosFreq,snr);

% Generate data realization
noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
dataVec = noiseVec+qcData;

% Input parameters for CRCBQCHRPPSO
inParams = struct('dataX', timeVec,...
                  'dataY', dataVec,...
                  'dataXSq',timeVec.^2,...
                  'dataXCb',timeVec.^3,...
                  'psdVal', psdPosFreq,...
                  'sampFreq',sampFreq',...
                  'rmin',rmin,...
                  'rmax',rmax);
              
% CRCBQCHRPPSO runs PSO on the CRCBQCHRPFITFUNC fitness function. As an
% illustration of usage, we change one of the PSO parameters from its
% default value.
outStruct = glrtqcsigpso(inParams,struct('maxSteps',2000),nRuns);

%%
% Plots
figure;
hold on;
plot(timeVec,inParams.dataY,'.');
plot(timeVec,qcData);
for lpruns = 1:nRuns
    plot(timeVec,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
end
plot(timeVec,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
legend('Data','Signal',...
       ['Estimated signal: ',num2str(nRuns),' runs'],...
       'Estimated signal: Best run');
disp(['Estimated parameters: a1=',num2str(outStruct.bestQcCoefs(1)),...
                             '; a2=',num2str(outStruct.bestQcCoefs(2)),...
                             '; a3=',num2str(outStruct.bestQcCoefs(3))]);

% PSO performance figure
figure;
hold on
for i=1:nRuns
   % plot every runs output in different color 
   % the differences of pso runs can be seen if one zoom into the graph
   plot(timeVec,outStruct.allRunsOutput(i).estSig,'SeriesIndex', i); %Requires R2020a or later version
   title("PSO Performance graph")
end
hold off
