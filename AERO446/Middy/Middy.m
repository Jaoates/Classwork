%% Joshua Oates - mid term test - 446
clear all
close all
clc
addpath("C:\joshFunctionsMatlab\")
%% Problem 1
clear
mu_e = 398600;%km^3/s^2
r_e = 6378;%km
z_geo = 35786;
r = r_e+z_geo;

theta2 = asind(r_e*sind(100)/r)
FOR2 = theta2*2

%% Problem 2
clear


%% Problem 3
clear
s = 4;
p = 8;

v = 3.6;
vp = v*s

Ahc = 48;
Ahp = Ahc*p

%% Problem 4
clear 
V = 3.6;
C = 48;
DOD = .5;
eta = 1;
t = 1;

P = V*C*DOD*eta/t

%% Problem 5
clear
T = [200,100];
T = T+273.15;
ratio = T(2)^4/T(1)^4
%% Problem 6
clear
r = 5;
A = pi*r^2;
lam = [1.7,2.5];%cm
lam = lam/100;
G = 4*pi*A./lam.^2;
G = 10*log10(G)

theta = 70*lam/(r*2)

z_geo = 35786;
theta = deg2rad(theta);
arclength = z_geo*theta



%% Problem 7 a
clear 
c = 299792458;%m/s
h =  6.62607015e-34; %J*Hz^-1;
q = 1.60217663e-19;% coulombs

R = .64;
lam = .8;%micrometers
lam = lam*1e-6; %meters
QE = (R/lam)*(h*c/q)

%% Problem 7 b
clear 
c = 299792458;%m/s
h =  6.62607015e-34; %J*Hz^-1;
q = 1.60217663e-19;% coulombs

R = .46;
lam = 1.1;%micrometers
lam = lam*1e-6; %meters
QE = (R/lam)*(h*c/q)

%% Problem 8
clear
C = 200;%Mbit
C = C*1e6;
B = 30;%MHz
B = B*1e6;

S_N = 2^(C/B)-1

SNR = 10*log10(S_N)

%% Problem 9
clear
syms t
A2_A1 = 2*pi/int(sin(t),[0,pi])





