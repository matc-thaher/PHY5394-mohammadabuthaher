function [fitVal,varargout] = crcbpsotestfunc_rosenbrock(xVec,params)
%A benchmark test function for CRCBPSO
%F = ATCCRCBPSOTESTFUNC(X,P)
%Compute the Rosenbrock fitness function for each row of X.  The fitness
%values are returned in F. X is standardized, that is 0<=X(i,j)<=1. P has
%two arrays P.rmin and P.rmax that are used to convert X(i,j) internally to
%actual coordinate values before computing fitness: X(:,j) ->
%X(:,j)*(rmax(j)-rmin(j))+rmin(j). 
%
%For standardized coordinates, F = infty if a point X(i,:) falls
%outside the hypercube defined by 0<=X(i,j)<=1.
%
%[F,R] =  CRCBPSOTESTFUNC(X,P)
%returns the real coordinates in R. 
%
%[F,R,Xp] = CRCBPSOTESTFUNC(X,P)
%Returns the standardized coordinates in Xp. This option is to be used when
%there are special boundary conditions (such as wrapping of angular
%coordinates) that are better handled by the fitness function itself.
%==========================================================================

%rows: points
%columns: coordinates of a point
[nrows,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nrows,1);
validPts = ones(nrows,1);

%Check for out of bound coordinates and flag them
validPts = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;
%Convert valid points to actual locations
xVec(validPts,:) = s2rv(xVec(validPts,:),params);


for lpc = 1:(nrows-1)
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);
        xnext = xVec(lpc + 1,:);
        %fitVal(lpc) = sum(100 .* (x(2:end) - x(1:end-1).^2).^2 +(x(1:end-1) - 1).^2);
        fitVal(lpc) = sum(100 .* (xnext - x.^2).^2 +(x - 1).^2);
        %fitVal(lpc) =4*x(1)^2 - 2.1*x(1)^4 + (x(1)^6)/3 + x(1) * x(2)- 4*x(2)^2 + 4*x(2)^4; % it works when x is being calculated as xVec(lpc,:)     
        %fitVal(lpc) =4.*x.^2 - 2.1.*x.^4 + (x.^6)/3 + x .* xnext - 4.*xnext.^2 + 4.*xnext.^4; % it works when x is being calculated as xVec(lpc)
    end
end

%Return real coordinates if requested
if nargout > 1
    varargout{1}=xVec;
    if nargout > 2
        varargout{2} = r2sv(xVec,params);
    end
end






