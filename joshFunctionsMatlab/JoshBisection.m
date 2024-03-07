function [r, count, output] = JoshBisection(f,bracket,TOL)
% this function (algorithm sourced from canvas) will use the bisection method to solve roots given a
% braket using Bisection Method
% f must have a single root on the interval bracket or this function will
% not return the expected root bracket
% JoshBisection is kinda janky but it works some of the time *shrug*

% Check inputs
arguments
    f
    bracket (1,2) {mustBeNumeric,mustBeReal}
    TOL (1,1) {mustBeNumeric,mustBeReal} = .1
end
if ~isa(f,'function_handle')
    throw(MException("JoshBisection:invalidInput","f must be a function_handle"))
end
% set vars
a = bracket(1);
b = bracket(2);

if ((sign(f(a)) == sign(f(b))) || (a > b))
    error('JoshBisection:invalidInput','f must not have the same sign on both sides of the bracket') 
end

% create output array
output = [a, b, f((a+b)/2)];
% start count = 0
count = 0;

% while err estimate > tol

isComplex = 0;
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
    if ~isreal(f(a)) || ~isreal(f(b)) || ~isreal(f(c))
        isComplex = 1;
    end

    if (f(a)*f(c) < 0)
        b = c;
    else
        a = c;
    end
    % update output array by appending new output
    output = [output; [a b f((a+b)/2)]];
end
% after err<=tol , return r is in middle of last range
if isComplex
    warning("JoshBisection: the function handle f returns complex numbers for tested values. JoshBisection was written for use on real functions only, but will still return a result that may be useful.")
end

r = (a+b)/2;
end