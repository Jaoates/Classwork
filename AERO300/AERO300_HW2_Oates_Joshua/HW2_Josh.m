%% section 0 - clean and clear
clear all
close all
clc

%% section 1 - solutions to matlab problems on hw2 for AERO 300 HW 2

%{
Problem 2 (40 points)
1. (17 points) Consider the equation x2 + ln x = 3 .
(a) (5 points) Find an interval [a,b] of length 1 that brackets the solution to the above equation.
    a) [1,2]

(b) (5 points) Using your answer from part (a), how many steps of the Bisection Method are
required to correctly approximate the solution within 10 decimal places?

    b) x2 + ln(x) - 3 = 0

    we know from n> (P+log10(b-a))/log10(2)
    that it will not take more than 34 iterations

(c) (7 points) Check your answer to part (b) using the bisection.m function from class (and now
on Canvas ). Be sure to show your published code!
    c) Experiementaly, it takes 29 iterations.
%}

nEstimate = (10+log10(2-1))/log10(2);
f = @(x) x^2 + log(x) - 3;
[r, count] = bisection(f,1,2,10e-10);

%{


2. (23 points) The 2D state of stress in a test specimen under combined loading is given by
σ =
[√2 −1
−1 0.5]
ksi.
%}

s=[[2^.5,-1];[-1,.5]];


%{
(a) (4 points) Determine the characteristic equation of the stress state σ .
Hint: Recall that the characteristic equation of a matrix A is the polynomial p given by
p(λ) = det(A −λI) ,
where I is the identity matrix.
1
%}
% P=det([[L*2^.5,-1];[-1,L*.5]]);
% P = .5 * 2^.5 - L* 2^.5 - L * .5 + L^2

%{
(b) (16 points) Use Matlab®and the Bisection method to determine the principal stresses in the
specimen to within 6 decimal places. Remember to publish your code.
Hint: The principal stresses (σ1 and σ2) are simply the roots of the characteristic polynomial
of the stress state σ .
%}

% lambda = L given by .5 * 2^.5 - L* 2^.5 - L * .5 + L^2 = 0
P = @(L) .5 * 2^.5 - L* 2^.5 - L * .5 + L^2;
[s1, c1] = bisection(P,0,1,10e-10);
[s2, c2] = bisection(P,1,2,10e-10);


%{
(c) (3 points) Using the Tresca failure criterion, predict whether the specimen will fail if the critical
shear stress of the material has been experimentally determined to be τcrit = 1.17 ksi.
Hint: The Tresca failure criterion is governed by the maximum shear stress in the specimen
which (in this case) can be taken to be τmax = (1/2)|σ1 −σ2|.
%}
Tmax = (1/2)*abs(s1-s2); % is less than 1.17 so the speciemen will not fail

%{
Problem 3 (40 points)
1. (4 points each) Find (by hand) all fixed points of the following g(x) :
(a) f(x) = (8 + 2x) / (2 + x^2)  = x
(b) f(x) = x^5 
    a) x = 2
    b) x = 1 and x = 0

2. (4 points each) Express each equation below as a fixed-point iteration problem x = g(x) in two
different ways:
(a) 1 + x −x^3 * e^x = 0 
(b) x^2 + ln x = 3 

FPI of
    a) 1 + 2x −x^3 * e^x = x     ; -1  + x^3 * e^x = x
    b) (ln x -3)/x = x           ; e^(x^2 - 3) = x

3. (24 points) In this part, we seek to approximate the principal stresses in the test specimen of Problem
2 part 2 using the fixed-point iteration method. As before, p(λ) is the characteristic polynomial of
the stress state as found in Problem 2 part 2(a).
(a) (3 points) Show that setting p(λ) = 0 is equivalent to setting 
(2λ^2 + √2 −2)/(1 + 2√2) = λ

Demonstrated using fpi
%}
g = @(x) (2*x^2 + 2^.5 -2)/(1 + 2* 2^.5);
TOL = 0.5*10^(-5);
[x1, c3] = fpi(g, 1.5 , TOL);



%{
(b) (6 points) Taking λ0 = 0.5 as the initial guess, use the function fpi.m to determine one of
the principal stresses in the specimen to within 6 decimal places. Can you determine the other
principal srress (for instance using a different initial guess)?

TOL = 0.5*10^(-5);
%}
[x2, c4] = fpi(g, .5 , 10e-6);
%{ 
you cannot find the other root using a different guess



(c) (12 points) Express p(λ) = 0 as a different fixed-point iteration problem, g(λ) = λ, and use
this new form to find the other principal stress in the specimen to within the same accuracy as
before.
.5 * 2^.5 - L* 2^.5 - L * .5 + L^2
=> -L^2=.5 * 2^.5 - L* 2^.5 - L * .5
=> L = -(.5*2^.5)/L + 2^.5 +.5
%}
g = @(x)  -(.5*2^.5)/x + 2^.5 +.5;

[x3, c3] = fpi(g, 1.4 , 10e-6);

%{

(d) (3 points) When compared to the bisection method used in Problem 2 part 2, which method
do you prefer? Which one was more efficient?
    d) i prefer bisection, its much cleaner and the predictability is good
    fpi is more effecient from a computational perspective but spent much
    more of my time
%}