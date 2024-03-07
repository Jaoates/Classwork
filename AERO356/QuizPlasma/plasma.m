clear; close all; clc

%% 3
TEC = 5.48e16; %/m2
z = 600; %km
f = 3e9; %hz
m_e = 9.11e-31; %kg
e = -1.602e-19; %C
eps0 = 8.854e-12; 
%K2 = e^2 / (8*pi^2*eps0*m_e);
K2 = 40.31;
c =  299792458;
k_b = 1.38e-23;

dr = K2*(TEC/(f^2))

%% 4

eta_p = 8e11; 
f_p = 8.98*sqrt(eta_p)*1e-6

%% 5
T = 1593; %K
% eta_p = 3.876e17;
eta_p = 3e11
lam = 4*pi/3 * (eps0*k_b*T / (eta_p^(1/3) * e^2))^(3/2)

%% 6

r = 2.2;
eta_p = 5e11;
z = 920;
Re = 6378;
a = z + Re;
v = sqrt(398600/a)*1e3;
A = pi*r^2;

I_i = A*-e*eta_p*v*1e3

%% 7
v = 591e3;
eta = 2.8* (.01)^(-3);
T = 1e6; %K
m_p = 1.67e-27; %kg
P_s = .5*m_p * eta * v^2 * 1e9
