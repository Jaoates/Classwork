%% housecleaning
clear all
close all
clc
addpath("C:\joshFunctionsMatlab\")

%% Problem2c
s = 7.4599;
v = 0.5226;
vofsof600 = joshInterp(7.3724,.4344,7.546,.4742,s);
mofv500to600 = joshInterp(500e3,.5226,600e3,vofsof600);
tau = (-1/v)*mofv500to600.m
a = sqrt(v/tau)

%% Problem2c units
u = symunit;
v = u.m^3/u.kg;
P = u.Pa;
m = -v/P; % negative bc slope is negative in the 
tau = (-1/v)*m;
a = sqrt(v/tau);
a = unitConvert(a,'SI','Base');
a = simplify(a)

%% Problem3
clear all

u = symunit;
M = 3;
T = 250*u.K;
P = 101*u.kPa;
rho = 1.4077*u.kg/u.m^3;
gamma = 1.4;
R = 287*u.J/(u.kg*u.K);
a = sqrt(gamma*R*T);
a = simplify(unitConvert(a,'SI','Base'));
P0 = P*(1+((gamma-1)/2)*M^2)^(1/(gamma-1));
T0 = T*(1+((gamma-1)/2)*M^2);
V = M*a

%% Problem5
clear all

u = symunit;

V1 = 500*u.m/u.s;
P1 = 50*u.kPa;
T1 = 250*u.K;

gamma = 1.4;
Cv = 0.718*u.kJ/(u.kg*u.K);
R = 287*u.J/(u.kg*u.K);
a1 = sqrt(gamma*R*T1);
a1 = simplify(unitConvert(a1,'SI','Base'));

M1 = V1/a1;
M2 = sqrt( (1+((gamma-1)/2)*M1^2) / (gamma*M1^2-(gamma-1)/2) );
T0 = T1*(1+((gamma-1)/2)*M1^2);
T2 = T0/(1+((gamma-1)/2)*M2^2)
V2 = M2*sqrt(gamma*R*T2);
V2 = simplify(unitConvert(V2,'SI','Base'))

% specific volume
v1 = simplify(unitConvert(R*T1/P1,'SI','Base'))
h0 = simplify(unitConvert(Cv*T1+P1*v1+V1^2/2,'SI','Base'))
P2v2 = h0 - Cv*T2 - V2^2/2
P2v2 = simplify(unitConvert(P2v2,'SI','Base'))

%% Problem 5 second try
clear all

gamma = 1.4;
R = 287;

V1 = 500;
P1 = 50e3;
T1 = 250;

a1 = sqrt(gamma*R*T1)
M1 = V1/a1

x1=1.560;
x2=1.580;
x=1.5776;


P2P1 = joshInterp(x1,.2673e1,x2,.2746e1,x);
T2T1 = joshInterp(x1,.1361e1,x2,.1374e1,x);
P02P01 = joshInterp(x1,.9097,x2,.9026,x);
M2M1 = joshInterp(x1,.6809,x2,.6746,x);

P01P1 = (1+M1^2*(gamma-1)/2)^(gamma/(gamma-1));
T0T1 = (1+M1^2*(gamma-1)/2)

T0 = T0T1*T1
P01 = P01P1*P1
T2 = T2T1*T1
P2 = P2P1*P1
M2 = M2M1*M1

a2 = sqrt(gamma*R*T2)
V2 = M2*a2

P02P2 = (1+M2^2*(gamma-1)/2)^(gamma/(gamma-1));

P02 = P02P2*P2

