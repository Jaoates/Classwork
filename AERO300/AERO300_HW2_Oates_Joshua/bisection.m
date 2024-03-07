function [r, count, output] = bisection(f, a, b, TOL)
% Bisection Method
% Given initial interval [a, b] such that f(a)*f(b) < 0
%Check for valid inputs to bisection method
if f(a)*f(b) >= 0 || a>b
    error('Either f(a)f(b)<0 is NOT satisfied or a>b !!!') %ceases execution
end
output = [a, b, f((a+b)/2)];
count = 0;
while ((b - a)/2 > TOL)
    count = count + 1;
    c = (a + b)/2;
    if f(c) == 0
        break;
    end
    if (f(a)*f(c) < 0)
        b = c;
    else
        a = c;
    end
    output = [output; [a b f((a+b)/2)]];
end
    
r = (a+b)/2;