%% final exam - joshua oaets

close all
clear all
clc

%% q1
g = 9.8;
U = 10;
L = (5/1000)*g;
cl_cd=5.8;
rho = 1.2;
q = .5*rho*U^2;

Cdi = (L/q)/cl_cd;

AR = 3.5;
A = 52/(100^2);
span = sqrt(AR*A);
mu = 1.72e-5;
Re = rho*U*span/mu;
d = span/AR;
CF = d*(1.328/(Re^.5))

dratio = Cdi/(Cdi+CF)
fratio = CF/(Cdi+CF)

