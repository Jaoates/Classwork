%%clear all
close all
clc

addpath('C:\AERO303\compressibleFlowRelations\')


%% 51.5
A =1.25^2;
Tinf = 30+273.15;
T1 = 90+273.15;

DT = T1-Tinf

k = .02808
Pr = .7202
nu = 1.896e-5
eps = 1

perimeter = 1.25*4
Lc = A/perimeter


Tf = mean([Tinf,T1])

Beta = 1/Tf

% Ra = Pr*g*Beta*(DT)*Lc^3/(nu^2)
g = 9.81
Gr = g*Beta*DT*Lc^3/nu^2
Ra = Gr*Pr

Nus = .27*Ra^.25

% Nus = h*Lc/k
h=k*Nus/Lc

q_con = A*h*DT

sig = 5.67e-8

q_rad = sig*eps*A*T1^4

F = 1

q_r = F*A*sig*(T1^4-Tinf^4)
% q_r = eps*A*sig*T1^4


%% Problem 62.5
clear 

% wheres the normal shock
gam = 1.4;

P1 = 348
Pe = 101.3

Pe_P1 = Pe/P1;
Ae_At = 2.637;

isen1 = joshComp(Ae_At,1,"AB");
pRatioMax = 1/isen1.P0_P

isen2 = joshComp(Ae_At,1,"AA");
pRatioNSexit = 1/isen2.P0_P
Me = isen2.M

NS1 = joshComp(Me,2,"M")
pRatioAcrossNS = NS1.P2_P1

pRatioMin = pRatioAcrossNS*pRatioNSexit

Pe = P1/isen2.P0_P



P02_P01 = Pe_P1*isen2.P0_P*

% x = 1/(Pe_P1*Ae_At);

% Me = -1/(gam-1)+sqrt(1/(gam-1)^2+(2/(gam-1))*(2/(gam+1))^((gam+1)/(gam-1))*(x)^2);
% Me = sqrt(Me);
% isen1 = joshComp(Me,1,"M");
% 
% P0e_Pe = isen1.P0_P;
% clear isen1
% 
% P0e_P1 = P0e_Pe*Pe_P1; % Poe_Poi
% 
% ns = joshComp(1/P0e_P1,2,"P0");
% % M2 = ns.M1;
% M1 = ns.M2;
% clear ns
% isen2 = joshComp(M1,1,"M");
% locNS = isen2.A_Ax
% 
% %part two
% r1 = joshComp(Ae_At,1,"AA");
% M1 = r1.M;
% r2 = joshComp(M1,2,"M");
% P02_P01 = r2.P02_P01;
% M2 = r2.M2;
% r3 = joshComp(M2,1,"M");
% P0e_Pe = r3.P0_P;
% Pe_P01_lower = P02_P01/P0e_Pe;
% r4 = joshComp(2,1,"AAs");
% Pe_P01_upper = 1/r4.P0_P




