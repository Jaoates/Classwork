%% Joshua Oates HW 5 gas dynamics
clear all
close all
clc

addpath("C:\AERO303\compressibleFlowRelations\")
addpath("C:\AERO303\CDNVirginaTech\")
%% problem 1
syms P0 A T0 M
gam = 1.4;
R = 287;
P = P0/(1+((gam-1)/2)*M^2)^(gam/(gam-1));
T = T0/(1+((gam-1)/2)*M^2);
rho = P/(R*T);
a = sqrt(gam*R*T);
V = M*a;

mdot = vpa(simplify(subs(rho*A*V,M,1)))
% mdot = vpa(simplify(rho*A*V))

%% problem 2
clear
gam = 1.4;
R = 287;

M1 = 2;
T1 = 250;
P01 = 75;

r1 = joshComp(M1,1,'M');
P0_P1 = r1.P0_P
T0_T = r1.T0_T;
T0 = T0_T*T1;
P01 = P0_P1*P01;

r2 = joshComp(M1,5,'M')

M2 = 2.385
r3 = joshComp(M2,1,'M')


P0_P2 = r3.P0_P
T0_T2 = r3.T0_T
P2 = P01*P0_P1/P0_P2
T2 = T0/T0_T2


r4 = joshComp(M2,5,'M')
mu2 = r4.mu

%% problem 3
clear 
Pi = 100;
Ti = 270;
Ae_Ai = 3;
M1 = 2;

r1 = joshComp(M1,1,'M')
Ans_Ax = r1.A_Ax*2
r2 = joshComp(Ans_Ax,1,'AA')
% P1 = Pi/(r2.P0_P)


r3 = joshComp(r2.M,2,'M')

P0i = r1.P0_P*Pi
T0i = r1.T0_T*Ti

P02 = P0i*r3.P02_P01
P1 = P02/r3.P02_P1
P2 = P1*r3.P2_P1

r4 = joshComp(r3.M2,1,'M')
Ae_Ax = (3/2)*r4.A_Ax
r5 = joshComp(Ae_Ax,1,'AB')

Pe = P02/r5.P0_P




%% problem 4
clear

M2 = 1.75;
P01 = 5;

r1 = joshComp(M2,1,'M')
r2 = joshComp(r1.A_Ax,1,'AB')
Pe = P01/r2.P0_P



% part 2
r3 = joshComp(M2,2,'M')
P02 = P01*r3.P02_P01
Pi = P02/r3.P02_P1
P2 = Pi*r3.P2_P1

% part 3
r1 = joshComp(M2,1,'M')
Pe = P01/r1.P0_P




%% Problem 5
clear 

% wheres the normal shock
gam = 1.4;

Pe_P1 = .7;
Ae_At = 2;
x = 1/(Pe_P1*Ae_At);


Me = -1/(gam-1)+sqrt(1/(gam-1)^2+(2/(gam-1))*(2/(gam+1))^((gam+1)/(gam-1))*(x)^2);
Me = sqrt(Me);
isen1 = joshComp(Me,1,"M");

P0e_Pe = isen1.P0_P;
clear isen1

P0e_P1 = P0e_Pe*Pe_P1; % Poe_Poi

ns = joshComp(1/P0e_P1,2,"P0");
% M2 = ns.M1;
M1 = ns.M2;
clear ns
isen2 = joshComp(M1,1,"M");
locNS = isen2.A_Ax

%part two
r1 = joshComp(Ae_At,1,"AA");
M1 = r1.M;
r2 = joshComp(M1,2,"M");
P02_P01 = r2.P02_P01;
M2 = r2.M2;
r3 = joshComp(M2,1,"M");
P0e_Pe = r3.P0_P;
Pe_P01_lower = P02_P01/P0e_Pe;
r4 = joshComp(2,1,"AAs");
Pe_P01_upper = 1/r4.P0_P


