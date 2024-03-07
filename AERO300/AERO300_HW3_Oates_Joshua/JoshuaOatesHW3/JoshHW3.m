% this is my matlab HW 3 for Aero 299
% Joshua Oates

%% Section 0 - cleanup
clear all;
close all;
clc;

%% Section 1 - Use Newtons on function

% create p of lambda for use in this lab
f =@(l) l^2-(2^.5+.5)*l+.5*2^.5-1;

% set up other vars for newtons and secant methods
TOL = 10e-6;
fp = @(l) 2*l-(2^.5+.5);

x0 = .1;
x1 = .2;
[r1N,c1N] = JoshNewtons(f,fp,x0);
[r1S,c1S] = JoshSecant(f,x0,x1);

x0 = 2.6;
x1 = 2.7;
[r2N,c2N] = JoshNewtons(f,fp,x0);
[r2S,c2S] = JoshSecant(f,x0,x1);

format long
disp("Newtons found roots at")
disp(r1N)
disp(r2N)
disp("in ")
disp(c1N)
disp(c2N)
disp("iterations.")
disp("Secant found roots at")
disp(r1S)
disp(r2S)
disp("in ")
disp(c1S)
disp(c2S)
disp("iterations.")

disp("in both cases, the functions found the roots in the same number of iterations")
