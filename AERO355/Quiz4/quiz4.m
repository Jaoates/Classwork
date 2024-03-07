%% Quiz 4 - particulates
clear all
close all
clc

%% Q1
dpart = .14;
rhopart = 5.8;
v = 13.6;
k = .8; % of target material
dcrater = k*dpart^(1.056)*rhopart^(.579)*v^(2/3)

%% Q2
mpart = .77;
rhopart = 5.6;
v = 4;
k = .32;
% critical thickness
tcrit = k*rhopart^(1/6)*mpart^(.352)*v^(.875)

%% Q3
t = 4;
A = 1.4;
F = .35;
P = 1-exp(-F*A*t)

%% Q4
t = 3;
A = .8;
F = .36;
PNCP = exp(-F*A*t)*100 % as percent

