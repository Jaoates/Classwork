%% Aero 421  - HW1 - Joshua Oates
clear all
close all
clc
addpath("C:\joshFunctionsMatlab\")

%% 12.1
mu = 398600;
mu = mu*1e9;

R0 = 7000;
tx = pi/4;
ty = pi/4;
tz = pi/4;
J = diag([100 120 80]);

[Cx,Cy,Cz]=joshAxisRotation("radian");

Cbg = Cx(tx)*Cy(ty)*Cz(tz)
Rg = [0 0 R0]'
rb = Cbg*Rg

Tg = (3*mu/(norm(rb)^5))*joshCross(rb)*J*rb

%% 13.1
close all
Ia = 150
It = 100

wz0 = .5
wx0 = .1
wy0 = .02
wt = sqrt(wx0^2+wz0^2)
ohm = -1.5*wz0
phi = acos(wy0/wt)

wx = @(t) wt*sin(ohm*t+phi)
wy = @(t) wt*cos(ohm*t+phi)
wz = @(t) ones(1,length(t))*wz0

figure 
hold on
t = linspace(0,10)
plot(t,wx(t))
plot(t,wy(t))
plot(t,wz(t))
legend("wx","wy","wz")

hz = Ia*wz0
ht = It*wt

h = sqrt(hz^2+ht^2)
gam = asin(ht/h)
ohmp = h/It



