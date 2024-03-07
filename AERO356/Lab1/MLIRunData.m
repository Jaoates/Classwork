%% MLI run data
% this file is for calculations and data output from the MLI lab while the 
% MLI was in the vacuum chamber
clear all
close all
load("MLIrun.mat")
t = MLIrun(:,1);
LT = MLIrun(:,2);
UT = MLIrun(:,3);
clear MLIrun

LT = LT + 273.15;
UT = UT + 273.15;
t = t*60;

f = figure;
hold on
plot(t,LT)
plot(t,UT)
legend("Bottom Temp","Top Temp","Interpreter","latex","Location","best")
xlabel("Time [s]",Interpreter="latex")
ylabel("Temp [K]",Interpreter="latex")
exportgraphics(f,'MLITempTime.png','Resolution',300)


e_m = .76; %mylar
e_k = .72; %kapton
ek = e_m; % interior layer
ei = e_k; % inner layer
eo = e_m; % outer layer
%   layers:
%          k | m t m t m h m t m t | m
%          i | n   n   n   n   n   | o

n = 5;
e_eff = (ei^-1+eo^-1-1+2*n*ek^-1-n)^-1 % effective

m = 1.03;
TBI = LT(1);
TTI = UT(1);
TBF = LT(end);
TTF = UT(end);
Tv = ((TBF-TBI)-(TTF-TTI))/m

A = 2*1.5;% in^2
A = A*6.4516;%cm^2
A = A/(100*100);
sig = 5.6704e-8;

T = LT(1);
dQ = e_eff*sig*A*T^4



