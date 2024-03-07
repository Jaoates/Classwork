%% 6.12
clear all
dia1 = 12.7;%mm
L = 50.8;%mm
L = L/1000;
P = 50;%kN
P = P/1000;
dia2 = 12.67494;%mm

dia2 = dia2/1000;
dia1 = dia1/1000;

E = 490/.007;

A = pi*(dia1/2)^2;
sig = P/A;

eps_a = sig/E;
eps_t = (dia1-dia2)/dia1;

nu = eps_t/eps_a