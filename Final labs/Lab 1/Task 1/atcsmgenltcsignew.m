function sigVec = atcsmgenltcsignew(timeData,timeSteps,snr, freq, phase)
% Generate a linear transient chirp signal
% S = ATCSMGENLTCSIG(T,TS,SNR,FQ,P)
% Generates a linear transient chirp signal S. T is the vector of
% time stamps at which the samples of the signal are to be computed. TS is 
% the 'struct' vector of time limits [t_a, t-a + L], signal has zero value 
% out of this interval. SNR is the matched filtering signal-to-noise ratio
% of S, FQ is the 'struct' vector of frequencies [f_0, f_1], and P is the
% initial phase. FQ and TS components that parametrize the phase of the 
% signal: f_0 * (t - t_a) + f_1 * (t - t_a).^2 + (phase / (2 * pi).

%Mohammad Abu Thaher Chowdhury, April 2021

idxt = find(timeData>=timeSteps.sTsig & timeData <= timeSteps.fTsig);
sigVec = zeros(1, length(timeData));
phaseVec = 2 * pi * (freq.inFreq * (timeData(idxt) - timeSteps.sTsig) + freq.fnFreq * (timeData(idxt) - timeSteps.sTsig).^2) + phase;
trigSin = sin(phaseVec);
sigVec(idxt) = snr * trigSin/norm(trigSin);

