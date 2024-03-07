%% Joshua Oates - 331 Midterm 1 - 2/9
clear all
close all
clc
addpath('C:\joshFunctionsMatlab\')

%% problem 1
clear 
syms x y z theta L
Ux = x*cos(theta*(z/L))-y*sin(theta*(z/L))-x
Uy = x*sin(theta*(z/L))-y*cos(theta*(z/L))-y
Uz = 0
diff =diff(Uy,x)

(theta*x*cos((theta*z)/L))/L + (theta*y*sin((theta*z)/L))/L;

%% problem 2

syms a rho t nu r ohm theta z
sig_rr = (rho*ohm^2/8)*(3+nu)*(a^2-r^2)
sig_tt = (rho*ohm^2/8)*((3+nu)*a^2-(3*nu+1)*r^2)

dsig_rrdr = diff(sig_rr,r)

sig = [[sig_rr,0,0];[0,sig_tt,0];[0,0,0]]
n = [a;0;0]
Tr = sig*n

sig_y = 250e6
rho = 7850
nu = .3
a = .5 
t = .1

ohm_T = sqrt((8/rho)*sig_y/(a^2*nu+3*a^2*nu))
ohm_M = sqrt((8/rho)*sig_y/(a^2*nu+3*a^2))

v = pi*a^2*t
m = v*rho
I = m*(a^2/2)
T = .5*I*ohm_M^2