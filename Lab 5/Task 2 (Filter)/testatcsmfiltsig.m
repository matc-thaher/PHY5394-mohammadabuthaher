%% Plot the filter of Siganl
% signal parameters
A_01 = 10;
A_02 = 5;
A_03 = 2.5;
f_01 = 100;
f_02 = 200;
f_03 = 300;
phase_01 = 0;
phase_02 = pi/6;
phase_03 = pi/4;

% number of samples
nsamples = 2048;

% sample frequency
samplFreq = 1024;

% Maximum Frequency is equal to nyquist frequency
maxFreq = samplFreq/2;

% Sampling time
timeData = (0:(nsamples-1))/nsamples;

% Generation of Signal
sigvec = atcsmfiltsig([A_01,A_02,A_03], [f_01, f_02, f_03], [phase_01, phase_02, phase_03], timeData);

%% Low pass Filter of the signal
%Filter order
fN = 30;
 
%Design the digital filter
filtDes = fir1(fN, (1.5 * f_01)/samplFreq, 'low');
 
% Applying low pass filter on signal
filtSig_l = fftfilt(filtDes, sigvec);

%% High pass Filter of the signal
%Design the digital filter
filtDes = fir1(fN, (0.7 * f_03)/samplFreq, 'high');
 
% Applying high pass filter on signal
filtSig_h = fftfilt(filtDes, sigvec);

%% Band pass Filter of the signal
%Design the digital filter
filtDes = fir1(fN,[(0.7 * f_02)/samplFreq (0.7 * f_03)/samplFreq]);
 
% Applying band pass filter on signal
filtSig_b = fftfilt(filtDes, sigvec);

%% Plots of filters
figure;

% Plot of Low Pass Filter
subplot(3,1,1)
plot(timeData, sigvec);
hold on
plot(timeData, filtSig_l);
hold off
legend(["Original Signal (s(t))", "Filtered Signal (s(t))"], "Location", "bestoutside")
ylabel('s(t)')
xlabel('Time in sec')
title("Low pass filter of signal")

% Plot of High Pass Filter
subplot(3,1,2)
plot(timeData, sigvec);
hold on
plot(timeData, filtSig_h);
hold off
legend(["Original Signal (s(t))", "Filtered Signal (s(t))"], "Location", "bestoutside")
ylabel('s(t)')
xlabel('Time in sec')
title("High pass filter of signal")
 
% Plot of Band Pass Filter
subplot(3,1,3)
plot(timeData, sigvec);
hold on
plot(timeData, filtSig_b);
hold off
legend(["Original Signal (s(t))", "Filtered Signal (s(t))"], "Location", "bestoutside")
ylabel('s(t)')
xlabel('Time in sec')
title("Band pass filter of signal")