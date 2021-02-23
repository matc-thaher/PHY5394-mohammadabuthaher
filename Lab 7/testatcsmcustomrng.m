%% Plotting the distribution

% Parameters of Uniform Distribution
a = -2;
b = 1;
% Parameters of Normal Distribution
m = 1.5;
s = 2.0;

% Call the function for Uniform Distribution Calculation 
resultUnfmDist = customrand(a,b);
% Call the function for Normal Distribution Calculation 
resultNrmlDist = customrandn(m,s);

%% PDF Calculation
% PDF calculation of Uniform Distribution
xRangeU = a-1: 1/length(resultUnfmDist): b + 1;
pdfUniform = zeros(1, length(xRangeU));
idxt = (xRangeU>= a & xRangeU<= b);
pdfUniform(idxt) = 1/(b-a);

%% Histogram and Plots of Distribution
% Plotting the Uniform Distribution
figure;
k = histogram(resultUnfmDist, 'normalization', 'pdf', 'BinWidth', 0.2);

% Plotting pdf curve on Uniform Distribution graph
hold on
plot(xRangeU, pdfUniform, 'LineWidth', 2)
hold off
title('Uniform Distribution')
xlabel("Result(values using random variable X")
ylabel("p_x(x) = U(x:a,b)")
legend(["Normalized Uniform Distribution", "Uniform PDF Curve"], 'Location', 'northeast')

% Plotting the Normal Distribution
figure;
histNormal = histogram(resultNrmlDist, 'normalization', 'pdf');

% Making pdf from histogram
pd_n = makedist('Normal', m, s);
xRangeN = histNormal.BinLimits(1): 1/length(resultNrmlDist) :histNormal.BinLimits(2);
pdfNormal = pdf(pd_n, xRangeN);

% Plotting pdf curve on Normal Distribution graph
hold on
plot(xRangeN, pdfNormal, 'LineWidth', 2)
hold off
title('Normal Distribution')
ylabel("p_x(x) = N(x:mu, sigma)")
xlabel("Result(values using random variable x")
legend(["Normalized Normal Distribution", "Normal PDF Curve"], 'Location', 'northeast')
