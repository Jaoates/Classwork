clear all;
close all;
clc

addpath("C:\joshFunctionsMatlab\")

h = 10*1000;
[T,P,rho,mu] = joshStdAtm(h);

gamma = 1.4;
Rgas = 287.85; % J kg K
a = sqrt(gamma*Rgas*T);

u= 210;
M = u/a;

T0 = T*(1+(gamma-1)/2)*M^2;

