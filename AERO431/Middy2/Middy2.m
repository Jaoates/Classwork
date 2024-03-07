%% Middy 2 - Joshua Oates
clear all
close all
clc

%% Problem 1
P = 2e3; %N
r = 36;% mm
r = r/1000;%m
t = 2;%mm
t = t/1000;
Apiston = pi*r^2;
Press = P/Apiston; % Pa
sig1 = Press*r/t % hoop stress
sig2 = Press*r/(t*2) % long stress

%% Problem 2

% applied torques
T1 = 3e3;%Nm
T2 = 11e3;%Nm

% L1 = L2 = .6
L = .6;

% torques through length
% T1 and 2 add together
T1 = T1+T2;

G = 75e9;%Pa
r = 50;%mm
r = r/1000;%m

J = pi*r^4/2; % polar moment of inertia

% U = T^2*L/(2*G*J)
U = (T1^2+T2^2)*L/(2*G*J) % internal strain energy


%% Problem 3

syms Nab Nabx Naby Nbc Nbd Nbdx Nbdy Nbe Nde Nae Pe Pa

theta1 = atan(2/1.5);
theta2 = atan(2/3);

eqn = [
Nbe == Pe
Naby == Pa
Nbdy == -Nbe-Naby
Nae==Nde
Nabx==-Nae
Nabx==-Nbdx-Nbc

Nabx==cos(theta1)*Nab
Naby==sin(theta1)*Nab

Nbdx==cos(theta2)*Nbd
Nbdy==sin(theta2)*Nbd
]

vars = [Nab Nabx Naby Nbc Nbd Nbdx Nbdy Nbe Nde Nae];
sol = solve(eqn,vars);

sol.Nbc = Pa*(4.5/2)+Pe*(3/2)

solu = subs(subs(sol,Pe,0),Pa,1)
soll = subs(subs(sol,Pe,60),Pa,40)
soll = subs(subs(sol,Pe,51e3),Pa,51e3*(2/3))

s1 = sqrt(2^2+1.5^2)
s2 = sqrt(2^2+3^2)
L = [s1,3,s2,3,1.5]'
n = double([solu.Nab solu.Nbc solu.Nbd solu.Nde solu.Nae]')
N = double([soll.Nab soll.Nbc soll.Nbd soll.Nde soll.Nae]')

A = 400;%mm
A = A/(1e3*1e3)
E = 200e9

del = sum(n.*N.*L)/(A*E)
del = del*1000

%%
3.24*29/18
