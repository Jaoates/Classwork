%% HW 6 - Props - Joshua Oates
clear all
close all
clc

%% Problem 1
g = 9.8;
cs = 1800;
gam = 1.52;
a = .0018;
n = .07;
rhofuel = 920;
rhoox = 1170;
mdotox = 2844;
t = 120;
grainOD = 3.8;
Pc0 = 4.8e6;
OFi = 2;
At = 1.6;
L = 40;
epsilon = 7;

dpi = .59;
% db = .34;
z = 50; %km
N = 7;

dp =  2*((a*(2*n+1)*(mdotox/(pi*N))^n*t+(dpi/2)^(2*n+1))^(1/(2*n+1)))
OF = mdotox^(1-n)*dp^(2*n-1)/(4^n*rhofuel*pi^(1-n)*a*L)
Ap = pi*(dp/2)^2
rdot = a*(mdotox/Ap)^n;
Ab = L*pi*dp*N;
mdotf = rhofuel*Ab*rdot;
mdottotal = mdotf+mdotox
Pc = mdottotal*cs/At

syms Pce
eqn1 = epsilon == (2./(gam+1)).^(1./(gam-1)).*Pce.^(1./gam).*(((gam+1)./(gam-1)).*(1-Pce.^((1-gam)./gam))).^-.5; % definition expansion ratio lecture3;
Pce = vpasolve(eqn1(1),Pce) ;
Pce = double(Pce) ;
Pec = 1/Pce;

Cf = sqrt(... % Coeffecient of thrust from lecture7
    ((2*gam.^2)./(gam-1)).*...
    ( (2./(gam+1)).^((gam+1)./(gam-1)) ).*...
    (1-Pec.^((gam-1)./gam)) )+...
    Pec*epsilon;
F = Cf*At*Pc
Isp = F/(mdottotal*g)

%% Problem 2
clear
g = 9.81;
gam = 1.2;
R = 259;
Tc = 3000;
Pc = 200e3;
At = .5e-4;
epsilon = 15;
mprop = 10;

syms Pce
eqn1 = epsilon == (2./(gam+1)).^(1./(gam-1)).*Pce.^(1./gam).*(((gam+1)./(gam-1)).*(1-Pce.^((1-gam)./gam))).^-.5; % definition expansion ratio lecture3
Pce = vpasolve(eqn1(1),Pce) ;
Pce = double(Pce) ;
Pec = 1/Pce;

Cf = sqrt((2*gam.^2./(gam-1)) .* ((2./(gam+1)).^ ((gam+1)./(gam-1)) .* (1-Pec.^((gam-1)./gam)))) + Pec.*epsilon ;

F = Cf*At*Pc
cs = sqrt(gam.*R.*Tc)./(gam.*sqrt((2./(gam+1)).^((gam+1)./(gam-1)))); % C star from lecture7

mdot = Pc*At/cs;
Isp = F/(mdot*g)
It = Isp*mprop*g
%% Problem 4
clear
d = 15e-2;
% P = 3e3;
M = 2.17e-25;
I = 2.2;
mprop = 20;
theta = 20;
theta = theta/2;
eta = .85;
Vb = 800;
q = 1.60217663e-19;
% vi = sqrt(2*q*Vb/M);
F = I*sqrt(2*M*Vb/q);
Ft = cosd(theta);
F = F*Ft
mdot = I*M/q;
mdot = mdot/eta
g = 9.81;
Isp = (Ft*eta/g)*sqrt(2*q*Vb/M)
t = mprop/mdot
It = Isp*mprop*g

%% Problem 5

dslot = [5,7,10,14]./100;
w = .1.*dslot;
rle = .5.*w;
m = 9.1093837e-31;
q = 1.60217663e-19;
Te = 3; 

B = ( 1./rle ) .*sqrt(m*8*Te/(pi*q))

