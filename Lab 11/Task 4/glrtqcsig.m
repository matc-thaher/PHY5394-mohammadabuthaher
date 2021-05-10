function glrt = glrtqcsig(timeVec, A, sampFreq, dataVec, psdPosFreq, inVec)


%Generate the unit norm signal (i.e., template). Here, the value used for
%'A' does not matter because we are going to normalize the signal anyway.
%Note: the GLRT here is for the unknown amplitude case, that is all other
%signal parameters are known.
sigVec = crcbgenqcsig(timeVec,A,inVec);
%We do not need the normalization factor, just the  template vector
[templateVec,~] = normsig4psd(sigVec,sampFreq,psdPosFreq,1);
% Calculate inner product of data with template
llr = innerprodpsd(dataVec,templateVec,sampFreq,psdPosFreq);
%GLRT is its square
glrt = llr^2;