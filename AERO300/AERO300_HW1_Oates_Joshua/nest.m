% From Numerial Analysis, Timothy Sauer (modified to ignore the basis
% points)
% Program 0.1 Nested Multiplication
% Evaluates the polynomial from nested form using Horner's Method
%Input: d - degree of polynomial
%       c - array of d+1 coefficients (constant term first)
%       x - point at which to evaluate the polynomials
%Output: y - value of polynomial evaluated at x
function y = nest(d, c, x)
y = c(d + 1);
for i = d:-1:1
    y = y.*x + c(i);
end