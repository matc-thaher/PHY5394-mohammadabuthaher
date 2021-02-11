function sigVec = atcsmfiltsig(snr, freq, phase, timedata)
%% Generate the signal of sinusoids
% S = ATCSMGENFILTSIG(SNR,FQ, P, TS)
% Generates a cumulative signal of sinusoid S. SNR is the matched filtering 
% signal-to-noise ratio of S. FQ is the vector of frequencies [f_01, f_02,
% f_03]. P is also a vecotr with [phase_01, phase_02, phase_03] and TS is 
% the vector of time stamps at which the samples of the signal are to be 
% computed.
% FQ and TS components that parametrize the phase of the signal:
% f_0 * (t - t_a) + f_1 * (t - t_a).^2 + (phase / (2 * pi). 

%Mohammad Abu Thaher Chowdhury, February 2021

% generation of phase vector for first signal
phaseVec_1 = (freq(1) * timedata) + phase(1)/ (2 * pi);

% generation of phase vector for second signal
phaseVec_2 = (freq(2) * timedata) + phase(2) / (2 * pi);

% generation of phase vector for third signal
phaseVec_3 = (freq(3) * timedata) + phase(3)/ (2 * pi);

% Generation of cumulative signal
sigVecSt = (snr(1) * sin ( 2 * pi * phaseVec_1)) + (snr(2) * sin ( 2 * pi * phaseVec_2)) + (snr(3) * sin ( 2 * pi * phaseVec_3));
sigVec = sigVecSt/norm(sigVecSt);