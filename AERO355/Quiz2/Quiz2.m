%% Joshua Oates Quiz 2 AREO 355 space environments
clear all
close all
clc

%% Solar constant
r0 = 1; %AU
r = 5.8; %AU
E0 = 1366.1; % W/m^2

E = E0*(r0/r)^2

%% Sun Spot Count
g = 4;
s = 25;
k = 1;
R = k*(10*g+s)



