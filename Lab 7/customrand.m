function resultUD = customrand(a,b,un)
% Generate a uniform distribution values for a certain range
% UD = CUSTOMRAND(A,B,UN)
% Generates a uniform distrbution values UD. A is the lower or starting
% limit of range and B is the upper or final limit of range. UN is the input
% parameter for rand function.

% Mohammad Abu Thaher Chowdhury, February 2021

x = rand(un(1), un(2));
resultUD = (b-a) .* x + a;


