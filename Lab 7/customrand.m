function resultUD = customrand(a,b,j,k)
% Generate a uniform distribution values for a certain range
% UD = CUSTOMRAND(A,B)
% Generates a uniform distrbution values UD. A is the lower or starting
% limit of range and B is the upper or final limit of range.

% Mohammad Abu Thaher Chowdhury, February 2021

x = rand(j,k);
resultUD = (b-a) .* x + a;


