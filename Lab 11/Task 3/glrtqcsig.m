function glrt = glrtqcsig( timeVec, A, sampFreq, dataVec, psdPosFreq, inVec)
% Generate a glrt
% GLRT = GLRTQCSIG(T, A, F, D, PSD, IN)
% Generate a glrt value for the data where GLRT is the result of it. T is a 
% time vector, A is amplitude, F is the sample Frequency, D is the data 
% Vector, PSD is the vector of Power spectral density values, and IN is the 
% input parameters [a1, a2, a3]. Here, the value used for
%'A' does not matter because we are going to normalize the signal anyway.

% Mohammad Abu Thaher Chowdhury, April 2021 
sigVec = crcbgenqcsig(timeVec,A,inVec);
%We do not need the normalization factor, just the  template vector
[templateVec,~] = normsig4psd(sigVec,sampFreq,psdPosFreq,1);
% Calculate inner product of data with template
llr = innerprodpsd(dataVec,templateVec,sampFreq,psdPosFreq);
%GLRT is its square
glrt = llr^2;
%disp(glrt)