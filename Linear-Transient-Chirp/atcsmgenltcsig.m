function sigVec = atcsmgenltcsig(timeData,timeSteps,snr,freq, phase)
% Generate a linear transient chirp signal
% S = ATCSMGENLTCSIG(T,TS,SNR,FQ,P)
% Generates a linear transient chirp signal S. T is the vector of
% time stamps at which the samples of the signal are to be computed. TS is 
% the vector of time limits [t_a, t-a + L], signal has zero value out of 
% this interval. SNR is the matched filtering signal-to-noise ratio of S,
% FQ is the vector of frequencies [f_0, f_1], and P is the initial phase.
% FQ and TS components that parametrize the phase of the signal:
% f_0 * (t - t_a) + f_1 * (t - t_a).^2 + (phase / (2 * pi). 

%Mohammad Abu Thaher Chowdhury, January 2021

idxt = find(timeData>=timeSteps(1) & timeData <= timeSteps(2));
sigVec = zeros(1, length(timeData));
phaseVec = freq(1) * (timeData(idxt) - timeSteps(1)) + freq(2) * (timeData(idxt) - timeSteps(1)).^2 + (phase / (2 * pi));
trigSin = sin(2 * pi * phaseVec);
sigVec(idxt) = snr * trigSin/norm(trigSin);

