%% Plotting the linear transient chirp signal
% Parameters and calculation for generating the signal
% Signal parameters
A = 5;
t_a = 3.0;
L = 4.0;
f_0 = 2;
f_1 = 5;
phase = 0;

% Maximum instantaneous frequency after t_a + L sec is
maxinstFreq = f_0 + 2 * f_1 * L;

% Sample Frequency with 5 times of the Maximum Frequency
samplFreq1 = 5 * maxinstFreq;
samplIntrvl = 1/samplFreq1;
% Sample Frequency with 1/2 times of the Maximum Frequency
samplFreq2 = 0.5 * maxinstFreq;
samplIntrv2 = 1/samplFreq2;

% Time samples for sample interval 1
timedata1 = 0:samplIntrvl:10;
% Time samples for sample interval 2
timedata2 = 0:samplIntrv2:10;

% Generate the signal for time data 1
sigVec1 = atcsmgenltcsignf(timedata1,[t_a, t_a + L], A,[f_0,f_1], phase);
% Generate the signal for time data 2
sigVec2 = atcsmgenltcsignf(timedata2,[t_a, t_a + L], A,[f_0,f_1], phase);

%% Plots of the signals
% Plotting for the first Nyquist Frequency
figure;
plot(timedata1,sigVec1, '* -')
xlabel("Time in sec")
ylabel("s(t)")
title(['Linear Transient Chirp Signal (Sampling frequency =',num2str(samplFreq1),' Hz)']);

%Plotting for the second Nyquist Frequency
figure;
plot(timedata2,sigVec2, '* -')
xlabel("Time in sec")
ylabel("s(t)")
title(['Linear Transient Chirp Signal (Sampling frequency =',num2str(samplFreq2),' Hz)']);


