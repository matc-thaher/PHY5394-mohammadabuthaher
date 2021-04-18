function sigVec = atcsmgenltcsignew(timeData,snr, sigParameters)
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

idxt = find(timeData>=sigParameters.sTsig & timeData <= sigParameters.fTsig);
sigVec = zeros(1, length(timeData));
phaseVec = 2 * pi * (sigParameters.inFreq * (timeData(idxt) - sigParameters.sTsig) + sigParameters.fnFreq * (timeData(idxt) - sigParameters.sTsig).^2) + sigParameters.phase;
trigSin = sin(phaseVec);
sigVec(idxt) = snr * trigSin/norm(trigSin);

