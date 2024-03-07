%% Props HW 3
clear all
close all
clc

addpath("C:\joshFunctionsMatlab\")
% addpath('C:\AERO303\compressibleFlowRelations\')
%% Problem 3
gam = 1.4
R = 287.05
MFP =@(M) sqrt(gam/R) * M * ( 1+(gam-1)*M^2/2 )^-( (gam+1)/(2*(gam-1)) );
fplot(MFP,[.1,2])
xlabel("Mach")
ylabel("MFP")
title("MFP over Mach for air")


%% Problem 4
u = symunit;

alt = 30000*u.ft;
alt = unitConvert(alt,u.m);
alt = double(separateUnits(alt));
[T,P,rho,mu] = joshStdAtm(alt);

M = .85
gam = 1.4

P_P0 = (1+((gam-1)/2)*M^2)^-(gam/(gam-1));
T_T0 = (1+((gam-1)/2)*M^2)^-1;
rho_rho0 = (1+((gam-1)/2)*M^2)^-(1/(gam-1));

% Tref = 518.67;
% Pref = 14.696;
% mdot = 102.56;
% mdotc = mdot*(sqrt(T/Tref))/(P/Pref)


P = P/P_P0
T = T/T_T0


syms gam M A R P_tot T_tot mdot
% gamma 
% Mach 
% AreaCompressor face
% P total
% T total
% mass flow rate
eqn = mdot ==( (A*P_tot)/sqrt(T_tot) ) * sqrt(gam/R) * M * ( 1+(gam-1)*M^2/2 )^-( (gam+1)/(2*(gam-1)) );


eqn = subs(eqn,mdot,102.56*u.lbm/u.s); % sub in mdot
eqn = subs(eqn,P_tot,P*u.Pa); % sub in P_tot
T_tot_in = T*u.K; % init T_tot_in
T_tot_in = unitConvert(T_tot_in,u.Rankine,"Temperature", "Absolute");
eqn = subs(eqn,T_tot,T_tot_in); % sub in T_tot
eqn = subs(eqn,gam,1.4); % sub in gam
R_in = 287 * u.J/(u.kg*u.K); % init R
R_in = unitConvert(R_in,"US"); % convert R 
R_in = unitConvert(R_in,u.Rankine,"Temperature", "Absolute"); % convert R
eqn = subs(eqn,R,R_in); % sub in R
eqn = subs(eqn,M,.5); % sub in M
A = simplify(solve(eqn,A)) % solve and simplify


syms gam M A R P_tot T_tot mdot
% gamma 
% Mach 
% AreaCompressor face
% P total
% T total
% mass flow rate
eqn = mdot ==( (A*P_tot)/sqrt(T_tot) ) * sqrt(gam/R) * M * ( 1+(gam-1)*M^2/2 )^-( (gam+1)/(2*(gam-1)) );

alt = 20000*u.ft;
alt = unitConvert(alt,u.m);
alt = double(separateUnits(alt));
[T,P,rho,mu] = joshStdAtm(alt);

M = .85
gam = 1.4

P_P0 = (1+((gam-1)/2)*M^2)^-(gam/(gam-1));
T_T0 = (1+((gam-1)/2)*M^2)^-1;
rho_rho0 = (1+((gam-1)/2)*M^2)^-(1/(gam-1));

P = P/P_P0
T = T/T_T0

syms gam M A R P_tot T_tot mdot
% gamma 
% Mach 
% AreaCompressor face
% P total
% T total
% mass flow rate
eqn = mdot ==( (A*P_tot)/sqrt(T_tot) ) * sqrt(gam/R) * M * ( 1+(gam-1)*M^2/2 )^-( (gam+1)/(2*(gam-1)) );


eqn = subs(eqn,mdot,102.56*u.lbm/u.s); % sub in mdot
eqn = subs(eqn,P_tot,P*u.Pa); % sub in P_tot
T_tot_in = T*u.K; % init T_tot_in
T_tot_in = unitConvert(T_tot_in,u.Rankine,"Temperature", "Absolute");
eqn = subs(eqn,T_tot,T_tot_in); % sub in T_tot
eqn = subs(eqn,gam,1.4); % sub in gam
R_in = 287 * u.J/(u.kg*u.K); % init R
R_in = unitConvert(R_in,"US"); % convert R 
R_in = unitConvert(R_in,u.Rankine,"Temperature", "Absolute"); % convert R
eqn = subs(eqn,R,R_in); % sub in R
eqn = subs(eqn,M,.5); % sub in M
A = simplify(solve(eqn,A)) % solve and simplify

%% Problem 5
% setup
clear
u = symunit;

% station 0 - infinity
Pt0 = 101.3*u.kPa % given
P0 = Pt0
Tt0 = 288*u.K % given
mdot0 = 77*u.kg/u.s

% station1 - opening of inlet
Pt1 = Pt0
Tt1 = Tt0
mdot1 = mdot0

% station2 - compressor face
Pt2 = Pt1
Tt2 = Tt1
mdot2 = mdot1

% station3 - after compressor, assume isentropic compression
PR = 13.5; % given
syms gam P_P0 T_T0 rho_rho0
assume(T_T0,"real")
eqn1 = P_P0 == T_T0^(gam/(gam-1));
eqn1 = subs(eqn1,P_P0,PR);
eqn1 = subs(eqn1,gam,1.4);
T_T0 = solve(eqn1,T_T0);

Tt3 = T_T0*Tt2
Pt3 = Pt2*PR
mdot3 = mdot2

% station4 - after combuster, assume no pressure loss
TET = 1260*u.K; % given
Tt4 = TET
Pt4 = Pt3
mdotfuel = 1.3*u.kg/u.s;
mdot4 = mdot3+mdotfuel

% station5 - after turbine, assume isentropic expansion
syms Tt5 Cph Cpc
Cph = 1150*u.J/(u.kg*u.K);
Cpc = 1005*u.J/(u.kg*u.K);
eqn = mdot4*Cph*(Tt4-Tt5) == mdot3*Cpc*(Tt3-Tt2); % governing eq
% eqn = subs(eqn,Cph,1150*u.J/(u.kg*u.K)); % sub in Cph
% eqn = subs(eqn,Cpc,1005*u.J/(u.kg*u.K)); % sub in Cpc
Tt5 = solve(eqn,Tt5)

syms Pt5 
eqn = P_P0 == T_T0^(gam/(gam-1)); % this can be used to find the total conditions across a turbine?
eqn = subs(eqn,P_P0,Pt5/Pt4);
eqn = subs(eqn,T_T0,Tt5/Tt4);
eqn = subs(eqn,gam,1.33); % sub in gamma hot
Pt5 = solve(eqn,Pt5)
mdot5 = mdot4

% station6 (after afterburner)- assume ideal combustion
mdotfuel = 1.9*u.kg/u.s;
LHV = 43.1*u.MJ/u.kg; % lower heating value or net calorific value
Qdotprime = mdotfuel*LHV; % ideal heat from fuel
mdot6 = mdot5+mdotfuel
Tt6 = simplify(Tt5+(Qdotprime/mdot6)/Cph)
Pt6 = Pt5

% converging nozzle - best possible is M = 1

mdot9 = mdot6
gam = 1.33;
Pt9 = Pt6
Tt9 = Tt6
M = 1; %ideal
P_P0 = (1+((gam-1)/2)*M^2)^-(gam/(gam-1));
T_T0 = (1+((gam-1)/2)*M^2)^-1;
P9 = P_P0*Pt9
T9 = T_T0*Tt9

R = 287*u.J/(u.kg*u.K);
V9 = M*sqrt(gam*R*T9)

syms Ae
eqn = mdot9 ==( (Ae*Pt9)/sqrt(Tt9) ) * sqrt(gam/R) * M * ( 1+(gam-1)*M^2/2 )^-( (gam+1)/(2*(gam-1)) );
Ae = simplify(solve(eqn,Ae));
Ae = unitConvert(Ae,u.m^2)
F = simplify(mdot9*V9 + (P9-P0)*Ae);
F = unitConvert(F,u.N)

TSFC = (3.2*u.kg/u.s)/F;
TSFC = unitConvert(TSFC,u.s/u.m)

% cd nozzle - best possible is p9 = p0
P9 = P0
P_P0 = P9/Pt9;

syms M    
assume(M,"positive")
eqn1 = P_P0 == (1+((gam-1)/2)*M^2)^-(gam/(gam-1));
M = solve(eqn1,M)
T_T0 = (1+((gam-1)/2)*M^2)^-1;
T9 = Tt9*T_T0
V9 = M*sqrt(gam*R*T9);
V9 = expand(V9);
V9 = unitConvert(V9,"SI");
V9 = simplify(V9)
syms Ae
eqn = mdot9 ==( (Ae*Pt9)/sqrt(Tt9) ) * sqrt(gam/R) * M * ( 1+(gam-1)*M^2/2 )^-( (gam+1)/(2*(gam-1)) );
Ae = (solve(eqn,Ae));

Ae = expand(Ae);
Ae = unitConvert(Ae,"SI");
% Ae = simplify(Ae);

Ae = unitConvert(Ae,u.m^2)
F = simplify(mdot9*V9 + (P9-P0)*Ae);
F = unitConvert(F,"SI","derived");
F = simplify(F);
F =  unitConvert(F,u.N)

TSFC = (3.2*u.kg/u.s)/F;
TSFC = unitConvert(TSFC,u.s/u.m)

%% mdot for ideal compressible gas - isentropic relations / nozzle relation
syms gam M A R P_tot T_tot mdot
% gamma 
% Mach 
% Area Compressor face / Area nozzle
% P total
% T total
% mass flow rate
eqn = mdot ==( (A*P_tot)/sqrt(T_tot) ) * sqrt(gam/R) * M * ( 1+(gam-1)*M^2/2 )^-( (gam+1)/(2*(gam-1)) );

% https://www.grc.nasa.gov/www/k-12/rocket/mflchk.html

%%  Pressure, Temperature, rho ratios for ideal compressible gas - isentropic relations/nozzle relation
syms gam P_P0 T_T0 rho_rho0
% gamma
% P/P_tot
% T/T_tot
% rho/rho_tot
eqn1 = P_P0 == rho_rho0^gam;
eqn2 = P_P0 == T_T0^(gam/(gam-1));
eqn3 = rho_rho0^gam == T_T0^(gam/(gam-1));

% https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html

%% Turbine work equation, class notes
syms mdot4 mdot3 Tt2 Tt3 Tt4 Tt5 Cph Cpc
% mass flow of air and fuel 
% mass flow of air
% Tt2 Temperature of gas entering compressor
% Tt3 Temperature of gas entering combustor
% Tt4 Temperature of gas entering turbine
% Tt5 Temperature of gas exiting turbine
% Specific heat of the gas exiting the combustor
% specific heat of the gas entering the combustor
eqn = mdot4*Cph*(Tt4-Tt5) == mdot3*Cpc*(Tt3-Tt2);

%% heat added to flow in combustor or afterburner
syms Qprime mdotfuel LHV
% Theoretical max heat from fuel [energy/time]
% mass flow of fuel [mass/time]
% lower heating value [energy/mass]
eqn1 = Qdotprime == mdotfuel*LHV;

syms Qdot mdotair Cph Tt3 Tt4
% mass flow of air
% specific heat of gas after combustor
% Tt3 Temperature of gas entering combustor
% Tt4 Temperature of gas entering turbine
eqn2 = Qdot == (mdotair+mdotfuel)*Cph*(Tt4-Tt3);

syms etaBurner 
% burner effeciency
eqn3 = etaBurner == Qdot/Qdotprime;

%% Mass Flow Parameter
syms gam M R MFP
% gamma 
% Mach 
% Area Compressor face / Area nozzle
% P total
% T total
% mass flow rate
eqn = MFP == sqrt(gam/R) * M * ( 1+(gam-1)*M^2/2 )^-( (gam+1)/(2*(gam-1)) );

%% Pressure Temperature and rho ratio and Mach - isentropic relations/nozzle relation
syms M gam P_P0 T_T0 rho_rho0
eqn1 = P_P0 == (1+((gam-1)/2)*M^2)^-(gam/(gam-1));
eqn2 = T_T0 == (1+((gam-1)/2)*M^2)^-1;
eqn3 = rho_rho0 == (1+((gam-1)/2)*M^2)^-(1/(gam-1));

% %% CPR - gives P/P0 for M=1. is useful for solving nozzles
% syms gam P_P0
% eqn = P_P0 == (2/(gam+1))^(gam/(gam-1))

%% Thrust
syms F mdot9 V9 P9 P0 Ae
% thrust
% mdot at exit
% velo at exit
% pressure at exit
% pressure at infinity distance
% Area of the exit
eqn = F == mdot9*V9 + (P9-P0)*Ae;


