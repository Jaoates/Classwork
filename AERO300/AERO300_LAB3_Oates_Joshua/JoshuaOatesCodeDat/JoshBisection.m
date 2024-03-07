function [r, count, output] = JoshBisection(f,braket,TOL)
% this function (algorithm sourced from canvas) will use the bisection method to solve roots given a
% braket using Bisection Method

% Check inputs
arguments
    f
    braket (1,2) {mustBeNumeric,mustBeReal}
    TOL (1,1) {mustBeNumeric,mustBeReal}
end

% set vars
a = braket(1);
b = braket(2);

if ((sign(f(a)) == sign(f(b))) || (a > b))
    error('invalid input') 
end

% create output array
output = [a, b, f((a+b)/2)];
% start count = 0
count = 0;

% while err estimate > tol
while ((b - a)/2 > TOL)
    % increment count
    count = count + 1;
    % create guess in middle of range
    c = (a + b)/2;
    % check if guess is root
    if f(c) == 0
        break;
    end
    % check if root is on left or right of guess
    if (f(a)*f(c) < 0)
        b = c;
    else
        a = c;
    end
    % update output array by appending new output
    output = [output; [a b f((a+b)/2)]];
end
% after err<=tol , return r is in middle of last range
r = (a+b)/2;
end