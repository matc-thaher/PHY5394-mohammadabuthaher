%% Test harness for CRCBPSO 
addpath 'D:\UTRGV_Spring_2021\Statistical_Methods\SDMBIGDAT19\CODES'
% The fitness function called is CRCBPSOTESTFUNC. 
ffparams = struct('rmin',-100,...
                     'rmax',100 ...
                  );
% Fitness function handle.
fitFuncHandle = @(x) crcbpsotestfunc_rosenbrock(x,ffparams);
%%
% Call PSO.
rng('default')
psoOut = crcbpso(fitFuncHandle,30);

%% Estimated parameters
% Best standardized and real coordinates found.
stdCoord = psoOut.bestLocation;
[~,realCoord] = fitFuncHandle(stdCoord);
disp(['Best location:',num2str(realCoord)]);
disp(['Best fitness:', num2str(psoOut.bestFitness)]);


%% Surface Plot 
% Rosenbrock Function
%[X,Y]=meshgrid(linspace(0,1,30));
%rsnbk=100*(Y-X.^2).^2+(1-X).^2;
%figure;
%surf(X,Y,rsnbk)

% Six hump camel function
%[X,Y]=meshgrid(linspace(-3,3,30));
%cam =4.*X.^2 - 2.1.*X.^4 + (X.^6)./3 + X .* Y - 4.*Y.^2 + 4.*Y.^4;
%figure;
%surf(X,Y,cam)
