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