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
fN = 50;
 
%Design the digital filter
filtDes = fir1(fN, (1.5 * f_01)/samplFreq, 'low');
 
% Applying low pass filter on signal
filtSig_l = fftfilt(filtDes, sigvec);

%% High pass Filter of the signal
%Design the digital filter
%filtDes = fir1(fN, (0.7 * f_03)/samplFreq, 'high');
filtDes = fir1(fN, (0.8 * f_03)/samplFreq, 'high');
 
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

%% Plot the Periodogram
%--------------
%Length of data 
dataLen = timeData(end)-timeData(1);
%DFT sample corresponding to Nyquist frequency
kNyq = floor(nsamples/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/dataLen);
% FFT of signal
fftSig = fft(sigvec);
fftSig_l = fft(filtSig_l);
fftSig_h = fft(filtSig_h);
fftSig_b = fft(filtSig_b);
% Discard negative Frequencies
fftSig = fftSig(1:kNyq);
fftSig_l = fftSig_l(1:kNyq);
fftSig_h = fftSig_h(1:kNyq);
fftSig_b = fftSig_b(1:kNyq);

%Plot periodogram
figure;

%plot of input signal
subplot(4,1,1)
plot(posFreq,abs(fftSig), 'm -')
xlabel("Frequency in Hz(only positive value)")
ylabel("DFT values of s(t)")
title("Discrete Fourier Transform of Input Signal")

%plot of periodogram of low pass filter
subplot(4,1,2)
plot(posFreq, abs(fftSig_l), 'g -')
xlabel("Frequency in Hz(only positive value)")
ylabel("DFT values of s(t)")
title("Discrete Fourier Transform of low pass filtered signal")

% Plot of periodogram of band pass filter
subplot(4,1,3)
plot(posFreq, abs(fftSig_h), 'r -')
xlabel("Frequency in Hz(only positive value)")
ylabel("DFT values of s(t)")
title("Discrete Fourier Transform of high pass filtered signal")

% Plot of periodogram of band pass filter
subplot(4,1,4)
plot(posFreq, abs(fftSig_b), 'b -')
xlabel("Frequency in Hz(only positive value)")
ylabel("DFT values of s(t)")
title("Discrete Fourier Transform of band pass filtered signal")
