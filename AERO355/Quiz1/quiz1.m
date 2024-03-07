%% quiz1
clear all
close all
clc

%% mean free path
rmol = 6.5e-10; %m
area = pi*rmol^2;
ng = 9.8e14; %1/m^3
% joshMeanFree(area,ng)
mfp = 1/(4*ng*area)

%% q0 first
T = 398; % K
TML = .9/100; %.7 percent
m = 2.7;%g
dm = TML*m; %g
Ea = 41.4; % kJ/mol
Ea = Ea*1000; % J/mol
R = 8.314; % J/(mol K) 
t2 = 24; % hrs
t1 = 0;

% T = 398; % K
% TML = 2.4/100; % 2.4 percent
% m = 2.2;%g
% dm = TML*m; %g
% Ea = 49.3; % kJ/mol
% Ea = Ea*1000; % J/mol
% R = 8.314; % J/(mol K) 
% t2 = 24; % hrs
% t1 = 0;

q0 = dm/(2*exp(-Ea/(R*T))*(t2^.5-t1^.5))

%% q0 second question
t2 = 24*4; % hrs
t1 = 0;
T = 207 ;
Ea = 31 ; % kJ/mol
Ea = Ea*1000; % J/mol
q0 = 22385020 ;
R = 8.314 ;
syms dm
eqn = q0 == dm/(2*exp(-Ea/(R*T))*(t2^.5-t1^.5));
dm = solve(eqn,dm);
dm = double(dm)

%% view factor
theta = 35;
phi = 52 ;
A1 = .03;
A2 =  .03;
r=2.2;
F = cosd(theta)*cosd(phi)*A2/(pi*r^2);
F100 = F*100



%% functions
% function [mfp,V,f,Kn] = joshMeanFree(area,ng,m,T,L)
%     % meanfree path
%     % Knudsen Number
%     % particle enrgy(optional)
%     % collison freq(optional)
%     % chracteristic length (optional)
%     % from collisional area of a particle, Temp, particle mass, number of particles in a
%     % m^3
%     arguments
%         area
%         ng
%         m = nan
%         T = nan
%         L = nan
%     end
%     K = 1.38e-23; % J/K
%     mfp = 1/(4*ng*area);
% 
%     if isnan(L) || isnan(T) || isnan(m)
%         Kn = nan;
%         V = nan;
%     else
%         V = sqrt(8*K*T/(pi*m));
%         Kn = mfp/L;
%         f = 4*V*ng*area;
%     end
% 
% end