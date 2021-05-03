%% Plot the linear transient chirp signal
% Adding path to the folder of function
addpath '../../Lab 1/Task 1'

% Signal parameters
%t_a = 0.2;
%L = 0.6;
%f_0 = 50;
%f_1 = 70;
%phase = 0;

% parameters
nSamples = 2048;
samplFreq = 2048;
% Time samples
timedata = (0:nSamples-1)/samplFreq;

% struct of the vectors
%sigTime = struct('sTsig', t_a, 'fTsig', t_a+L);
%sigFreq = struct('inFreq', f_0, 'fnFreq', f_1);
%sigParameters = struct('sTsig', t_a, 'fTsig', t_a + L, 'inFreq', f_0, 'fnFreq', f_1, 'phase', 0);
sigParameters = struct('sTsig', 0.2, 'fTsig', 0.8, 'inFreq', 50, 'fnFreq', 70, 'phase', 0);
snr = struct('snr1', 10, 'snr2', 12, 'snr3', 15);

% Generate the signal using function handle
H = @(snr) atcsmgenltcsignew(timedata,snr,sigParameters);

%% Plot the time series of the signal 
figure;
subplot(3,1,1)
plot(timedata,H(snr.snr1), '* -')
xlabel("Time in sec")
ylabel("signal")
title("Linear Transient Chirp Signal")
%subtitle("Linear Transient Chirp Signal with snr = 10")
subplot(3,1,2)
plot(timedata,H(snr.snr2), 'o -')
xlabel("Time in sec")
ylabel("signal")
%subtitle("Linear Transient Chirp Signal with snr = 12")
subplot(3,1,3)
plot(timedata,H(snr.snr3), '+ -')
xlabel("Time in sec")
ylabel("signal")
%subtitle("Linear Transient Chirp Signal with snr = 15")

