function sigVec = atcsmgenltcsignf(timeData,timeSteps,snr,freq, phase)
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

phaseVec = zeros(1,length(timeData));
trigSin = zeros(1,length(timeData));
sigVec = zeros(1,length(timeData));
for x = 1: length(timeData)
    if timeData(x) >= timeSteps(1) && timeData(x) <= timeSteps(2)
        phaseVec(x) = freq(1) * (timeData(x) - timeSteps(1)) + freq(2) * (timeData(x) - timeSteps(1)).^2 + (phase / (2 * pi));
        trigSin(x) = sin(2 * pi .* phaseVec(x));
        sigVec(x) = snr .* trigSin(x)/(norm(trigSin(x)));
    end
end

