Cp%% HW3 - Aero 303 - gas dynamics - Joshua Oates
clear all
close all
clc
addpath('C:\joshFunctionsMatlab\')
addpath('C:\AERO303\compressibleFlowRelations\')

%% Problem 1
clc
gam = 1.4;
R = .287;%kJ/(kg*K)
R = R*1000;%J/(kg*K)

cp = 1.005;%kJ/(kg*K)
cp = cp*1000;%J/(kg*K)

Pa_atm = 101325;

% state 1
P1 = 0.8;%atm
P1 = P1*Pa_atm;%Pa
M1 = 0.4;
T1 = 273;%K

% heat added
q = .25e6;%J/kg

% ideal gas
rho1 = P1/(R*T1);

% isentropic relations
T01 = T1*(1+M1^2*(gam-1)/2);
P01 = P1*(1+M1^2*(gam-1)/2)^(gam/(gam-1));
rho01 = rho1*(1+M1^2*(gam-1)/2)^(1/(gam-1));

% ratios at M1
P_Px = .1961e1;
T_Tx = .6151e0;
rho_rhox = .3188e1;
P0_P0x = .1157e1;
T0_T0x = .5290e0;

% from ratios
Px = P1/P_Px;
Tx = T1/T_Tx;
rhox = rho1/rho_rhox;
% T0x = T01/T0_T0x;
P0x = P01/P0_P0x;

T02 = q/cp + T01;
T02_Tx = (T02/T01)*(T0_T0x);

M2 = .9285;

% ratios at M2
P_Px = 1.0874;
T_Tx = 1.0195;
rho_rhox = 1.0666;
P0_P0x = 1.0025;
T0_T0x = .9961;

% from ratios
P2 = P_Px*Px;
T2 = T_Tx*Tx;
rho2 = rho_rhox*rhox;
P02 = P0_P0x*P0x;

P01 = P01/Pa_atm;
P1 = P1/Pa_atm;
Px = Px/Pa_atm;
P0x = P0x/Pa_atm;

M2
P2 = P2/Pa_atm
T2
rho2

T02 
P02 = P02/Pa_atm

%% Problem 2
clear
gam = 1.4;
R = .287;%kJ/(kg*K)
R = R*1000;%J/(kg*K)

cp = 1.005;%kJ/(kg*K)
cp = cp*1000;%J/(kg*K)

% state 1 givens
M1 =.2;
P1 = 100*1000;
T1 = 300;
q = 50e3;% J/(kg*K)

% state 1 solves
rho1 = P1/(R*T1);
T01 = T1*(1+M1^2*(gam-1)/2);
P01 = P1*(1+M1^2*(gam-1)/2)^(gam/(gam-1));
rho01 = rho1*(1+M1^2*(gam-1)/2)^(1/(gam-1));


% a call to compressible. input is M1, table is 3 (heat addition), input
% type is Mach #
out = compressible(M1,3,'M');
ratios.M = out(1);
ratios.P_Px = out(2);
ratios.T_Tx = out(3);
ratios.rho_rhox = out(4);
ratios.P0_P0x = out(5);
ratios.T0_T0x = out(6);
r1 = ratios;
clear out ratios

% from ratios star conditions
Px = P1/r1.P_Px;
Tx = T1/r1.T_Tx;
rhox = rho1/r1.rho_rhox;
T0x = T01/r1.T0_T0x;
P0x = P01/r1.P0_P0x;

% from heat transfer
T02 = q/cp + T01;
T02_T0x = (T02/T01)*(r1.T0_T0x);


% a call to compressible. input is T02_Tx, table is 3 (heat addition), input
% type is T0/T*
out = compressible(T02_T0x,3,'T0B');
ratios.M = out(1);
ratios.P_Px = out(2);
ratios.T_Tx = out(3);
ratios.rho_rhox = out(4);
ratios.P0_P0x = out(5);
ratios.T0_T0x = out(6);
r2 = ratios;
clear out ratios

P2 = r2.P_Px*Px;
T2 = r2.T_Tx*Tx;
rho2 = r2.rho_rhox*rhox;
P02 = r2.P0_P0x*P0x;
M2 = r2.M;

P01 = P01/1000;
P1 = P1/1000;

clc
M2
P2 = P2/1000
T2
rho2

T02 
P02 = P02/1000

%% Problem 3 
clear
clc
gam = 1.4;
R = .287;%kJ/(kg*K)
R = R*1000;%J/(kg*K)

cp = 1.005;%kJ/(kg*K)
cp = cp*1000;%J/(kg*K)

% state1 givens
M1=.6;
P1 = 150e3;%Pa
T1 = 300;%K
L = 45/100;%m
d = 3/100;%m
f = .005; % friction factor

fourFL_D = 4*f*L/d;

% state 1 solves
rho1 = P1/(R*T1);
T01 = T1*(1+M1^2*(gam-1)/2);
P01 = P1*(1+M1^2*(gam-1)/2)^(gam/(gam-1));
rho01 = rho1*(1+M1^2*(gam-1)/2)^(1/(gam-1));

% a call to compressible. input is M1, table is 4 (with friction), input
% type is Mach #
out = compressible(M1,4,'M');
ratios.M = out(1);
ratios.T_Tx = out(2);
ratios.P_Px = out(3);
ratios.rho_rhox = out(4);
ratios.P0_P0x = out(5);
ratios.fourFLx_D = out(6);
r1 = ratios;
clear out ratios

% from ratios star conditions
Px = P1/r1.P_Px;
Tx = T1/r1.T_Tx;
rhox = rho1/r1.rho_rhox;
P0x = P01/r1.P0_P0x;

% from gov eq
r2.fourFLx_D = r1.fourFLx_D - fourFL_D;

% a call to compressible. input is fourFLx_D2, table is 4 (with friction)
out = compressible(r2.fourFLx_D,4,'F');
ratios.M = out(1);
ratios.T_Tx = out(2);
ratios.P_Px = out(3);
ratios.rho_rhox = out(4);
ratios.P0_P0x = out(5);
ratios.fourFLx_D = out(6);
r2 = ratios;
clear out ratios

M2 = r2.M;
P2 = r2.P_Px*Px;
T2 = r2.T_Tx*Tx;

P01 = P01/1000;
P1 = P1/1000;

clc
M2
P2 = P2/1000
T2

%% Problem 4
clear
clc
gam = 1.4;
R = .287;%kJ/(kg*K)
R = R*1000;%J/(kg*K)

cp = 1.005;%kJ/(kg*K)
cp = cp*1000;%J/(kg*K)

% givens
V1 = 100; % m/s
T1 = 400; 
HV = 40e6;% J/kg

% state1 solve
a1 = sqrt(gam*R*T1);
M1 = V1/a1;
T01 = T1*(1+M1^2*(gam-1)/2);

% table call
out = compressible(M1,3,'M');
ratios.M = out(1);
ratios.P_Px = out(2);
ratios.T_Tx = out(3);
ratios.rho_rhox = out(4);
ratios.P0_P0x = out(5);
ratios.T0_T0x = out(6);
r1 = ratios;
clear out ratios

% from ratios star conditions
% Px = P1/r1.P_Px;
Tx = T1/r1.T_Tx;
% rhox = rho1/r1.rho_rhox;
T0x = T01/r1.T0_T0x;
% P0x = P01/r1.P0_P0x;

% table call
out = compressible(1,3,'M');
ratios.M = out(1);
ratios.P_Px = out(2);
ratios.T_Tx = out(3);
ratios.rho_rhox = out(4);
ratios.P0_P0x = out(5);
ratios.T0_T0x = out(6);
r2 = ratios;
clear out ratios

% from ratios
T02 = r2.T0_T0x*T0x;

% from governing eq
qmax = cp*(T02-T01);
fuel_airmax = cp*(T02-T01)/HV;

% if q was %10 higher
fuel_air=fuel_airmax*1.1;
q=qmax*1.1;

% clear M1 r1 r2 T02 T1 Tx V1

% new T0x based on new q
T0x = q/cp+T01;
T0_T0x=T01/T0x;

% table call
out = compressible(T0_T0x,3,'T0B');
ratios.M = out(1);
ratios.P_Px = out(2);
ratios.T_Tx = out(3);
ratios.rho_rhox = out(4);
ratios.P0_P0x = out(5);
ratios.T0_T0x = out(6);
r3 = ratios;
clear out ratios

% table call
% out = compressible(1,3,'M');
% ratios.M = out(1);
% ratios.P_Px = out(2);
% ratios.T_Tx = out(3);
% ratios.rho_rhox = out(4);
% ratios.P0_P0x = out(5);
% ratios.T0_T0x = out(6);
% r4 = ratios;
% clear out ratios

M3 = r3.M
% rho1 = P1/(R*T1);
T01_T3 = (1+M3^2*(gam-1)/2)
% P01 = P1*(1+M1^2*(gam-1)/2)^(gam/(gam-1));
rho01_rho3 = (1+M3^2*(gam-1)/2)^(1/(gam-1))
rho01_rho1 = (1+M1^2*(gam-1)/2)^(1/(gam-1))

T3 = T01/T01_T3
a3 = sqrt(gam*R*T3)
V3 = r3.M*a3
rho3_rho1 = rho01_rho1/rho01_rho3

m3_m1 = rho3_rho1*(V3/V1)
