% HW 4 - pressure vessels - Joshua Oates
%% Problem 1
clear all

P = 2000;%N
d = 47;%mm
d = d/1000;%m
r = d/2;
A = pi*r^2;%m^2
P = P/A;%Pa
t = 2;%mm
t = t/1000;%m

sig_hoop = P*r/t
sig_long = P*r/(2*t)


%% Problem 2
clear all

t = 8;%mm
t = t/1000;
d = 10;%mm
d = d/1000;
s = 50;%mm
s = s/1000;
P = 1.35;%Mpa
P = P*1e6;
r = .75;

sig_hoop1 = P*r/t
sig_hoop2 = P*r/(2*t);
sig_hoop2 = sig_hoop2*(50/40)
F = sig_hoop2*t*(40/1000)
A = pi*(d/2)^2;
shear = F/A


