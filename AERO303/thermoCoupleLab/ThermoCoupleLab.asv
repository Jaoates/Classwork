%% thermo couple lab
clear all
close all
clc

d1 = importdata("cooling.txt");
d2 = importdata("heating.txt");
d3 = importdata("solidrocket.txt");

dat1 = d1.data;
dat2 = d2.data;
dat3 = d3.data;

%%
f = figure;
hold on
plot(dat1(:,1),dat1(:,2))
plot(dat1(:,1),dat1(:,3))
plot(dat1(:,1),dat1(:,4))
plot(dat2(:,1),dat2(:,2))
plot(dat2(:,1),dat2(:,3))
plot(dat2(:,1),dat2(:,4))
ylabel("temp [deg F]","Interpreter","latex")
xlabel("time [s]","Interpreter","latex")

title("Temperature vs Time for 3 Materials Heating and Cooling","Interpreter","latex")
legend([["heating "+["Cu","Al","Steel"]],["cooling "+["Cu","Al","Steel"]]],"Location","best","Interpreter","latex")
f.Position = [f.Position]+200*[-1 -1 1 1];


%%
f= figure ;
hold on
plot(dat3(:,1),dat3(:,2))
plot(dat3(:,1),dat3(:,3))
plot(dat3(:,1),dat3(:,4))
plot(dat3(:,1),mean(dat3(:,2:4),2))
ylabel("temp [deg F]","Interpreter","latex")
xlabel("time [s]","Interpreter","latex")

title("Temperature vs Time for 3 Thermocouples During Estes firing","Interpreter","latex")
legend(["Temp 1","Temp 2","Temp 3","Temp Ave"],"Location","best","Interpreter","latex")
f.Position = [f.Position]+200*[-1 -1 1 1];

%% presentation
M = 1
T0 = 1916.76
P0 = 758.425
gam = 1.25
R = 287
Tin = 1273
Tout =  289.817 
mu = 48.445e-6
L = 0.0015875
% mu = joshInterp(273,17.456,293,18.24,T)

Ts = T0/(1+(gam-1)/2)
Ps = P0/(1+(gam-1)/2)^(gam/(gam-1))
rho = (Ps*1000)/(R*Ts)
u = M*sqrt(gam*R*T)

Re = 2*L*u*rho/mu

n = .4

Pr = .74
ratio = .00236
Nu = .023*Re^.8*Pr^n

kcg = .36
D = 2*L

h_bar = kcg*Nu/D

ka = .18
kb = .15

r1 = 1/16
r2 = .2825
r3 = .3450

Rt = 1/(2*pi*r1*L*h_bar)+log(r2/r1)/(2*pi*L*ka)+log(r3/r2)/(2*pi*L*kb)

DT = Tin - Tout

qburn = DT/Rt

qend = (1500-Tout)/(log(r2/r1)/(2*pi*L*ka)+log(r3/r2)/(2*pi*L*kb))
tburn = 1.6
ttot = 82
trat =  tburn/ttot
qave = qburn*(trat)+qend*(1-trat)

