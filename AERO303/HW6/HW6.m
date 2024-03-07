%% HW6 - A303 - Joshua Oates
clear all
close all
clc

u = symunit;
%% Problem 1
clear all
u = symunit;
L = 2*u.m;
d = .3*u.cm;

A = pi*d*L;

Tinf = (15+273.15)*u.K;
Ts = (152+273.15)*u.K;
% Tinf = unitConvert(Tinf,u.K);
% Ts = unitConvert(Ts,u.K);

Volts = 60*u.V;
I = 1.5*u.A;
q = Volts*I;
hc = q/(A*(Ts-Tinf));
hc = unitConvert(vpa(hc),u.W/(u.m^2*u.K))

%% Problem 2
clear all
u = symunit;

Tinf = (22+273.15)*u.K;
Twalls1 = (10+273.15)*u.K;
Twalls2 = (25+273.15)*u.K;
Tperson = (30+273.15)*u.K;

A = 1.4*u.m^2;
eps = .95;
sig = 5.6704e-8*u.W/(u.m^2*u.K^4);

q1 = eps*sig*A*(Tperson^4-Twalls1^4);
q2 = eps*sig*A*(Tperson^4-Twalls2^4);

q1 = vpa(unitConvert(q1,u.W))
q2 = vpa(unitConvert(q2,u.W))

%% Problem 3
clear all
L = .01;
T1 = 300+273.15
T2 = 200+273.15
DT = T1-T2
sig = 5.6704e-8
% a
k1 = .0219;
Rr = DT/(sig*(T1^4-T2^4));
Rk = L/k1;
R1 = Rr*Rk/(Rr+Rk);
q1 = DT/R1

% b
Rr = DT/(sig*(T1^4-T2^4));
q2 = DT/Rr


% c 
k2 = .026;
Rk = L/k2;
q3 = DT/Rk
% b
k3 = .00002;
Rk = L/k3;
q4 = DT/Rk




%% Problem 4
clear all
d1 = .05
d2 = .055
d3 = .055+.03
h1 = 60
h2 = 18
k1 = 80
k2 = .05

R1 = (1/(h1*pi*d1))
R2 = log(d2/d1)/(2*pi*k1)
R3 = log(d3/d2)/(2*pi*k2)
R4 = (1/(h2*pi*d3))



T1 = 320+273.15

R = (1/(h1*pi*d1))+(1/(h2*pi*d3))+log(d2/d1)/(2*pi*k1)+log(d3/d2)/(2*pi*k2)

DT = 320-5
q = DT/R

T2 = T1 - R1*q
T3 = T2 - R2*q
T4 = T3 - R3*q
DTp = T3-T2
DTin = T4-T3
