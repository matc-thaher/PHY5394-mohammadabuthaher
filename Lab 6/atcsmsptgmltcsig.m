function sigVec = atcsmsptgmltcsig(timeData,timeSteps,snr,freq, phase)
% Generate a linear transient chirp signal
% S = CRCBGENLSIG(X,SNR,C)
% Generates a linear transient chirp signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% three coefficients [a1, a2, a3] that parametrize the phase of the signal:
% a1*t+a2*t^2+a3*t^3. 

%Mohammad Abu Thaher Chowdhury, January 2021

idxt = find(timeData>=timeSteps(1) & timeData <= timeSteps(2));
sigVec = zeros(1, length(timeData));
phaseVec = freq(1) * (timeData(idxt) - timeSteps(1)) + freq(2) * (timeData(idxt) - timeSteps(1)).^2 + (phase / (2 * pi));
trigSin = sin(2 * pi * phaseVec);
sigVec(idxt) = snr * trigSin/norm(trigSin);
