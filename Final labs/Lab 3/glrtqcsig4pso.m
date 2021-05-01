function [fitVal,varargout] = glrtqcsig4pso(xVec, params)
%Fitness function for quadratic chirp regression
%F = GLRTQCSIG4PSO(X,P)
%Compute the fitness function (sum of squared residuals function after
%maximimzation over the amplitude parameter) for data containing the
%quadratic chirp signal. X.  The fitness values are returned in F. X is
%standardized, that is 0<=X(i,j)<=1. The fields P.rmin and P.rmax  are used
%to convert X(i,j) internally before computing the fitness:
%X(:,j) -> X(:,j)*(rmax(j)-rmin(j))+rmin(j).
%The fields P.dataY and P.dataX are used to transport the data and its
%time stamps. The fields P.dataXSq and P.dataXCb contain the timestamps
%squared and cubed respectively.
%
%[F,R] = GLRTQCSIG4PSO(X,P)
%returns the quadratic chirp coefficients corresponding to the rows of X in R. 
%
%[F,R,S] = GLRTQCSIG4PSO(X,P)
%Returns the quadratic chirp signals corresponding to the rows of X in S.
%==========================================================================

%rows: points
%columns: coordinates of a point
[nVecs,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nVecs,1);

%Check for out of bound coordinates and flag them
validPts = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;
for m = 1:3
    xVec(validPts,m) = xVec(validPts,m)*(params.rmax(m)-params.rmin(m))+params.rmin(m);
    %s2rv(xVec(validPts,:),params);
end
for lpc = 1:nVecs
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);
        fitVal(lpc) = ssrqc(x, params);
    end
end

%Return real coordinates if requested
if nargout > 1
    varargout{1}=xVec;
end

%log likelihood ratio after maximizing over amplitude parameter
function ssrVal = ssrqc(x, params)
%Generate normalized quadratic chirp
phaseVec = x(1).*params.dataX + x(2).*params.dataXSq + x(3).*params.dataXCb;
qc = sin(2.*pi.*phaseVec);

% Compute norm for colored data
[~,normFac] = normsig4psd(qc, params.sampFreq, params.psdVal, params.snr);
% Colored QC
qc =  qc./normFac;

% Compute fitness
ssrVal = -(innerprodpsd(params.dataY,qc,params.sampFreq,params.psdVal))^2;