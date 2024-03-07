%% [AERO452 PROJECT 1: GOES-18] Josh Oates, David Kumar
clear all;close all;clc;
load('TLE.mat')

mu = 398600; 
n = tle.MeanMotion;
n = deg2rad(n);
T = 2*pi/n; % sec
MeanAnomaly = tle.MeanAnomaly;

a = 42164; % km
ecc = 10^-9; % zero
incl = 10^-9; % zero
RAAN = tle.RightAscensionOfAscendingNode;
argp = tle.ArgumentOfPeriapsis;
nu = 0;


[r_vect, v_vect] = keplerian2ijk(a, ecc, incl, RAAN, argp, nu);

r_vect = [3207 5459 2714]; % example r and v vectors
v_vect = [ -6.532 0.7835 6.142];

state = [r_vect,v_vect];
t = [0 24*3600];
muearth = 398600;

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,statenew]=ode45(@pig,t,state,options,muearth);

figure
plot3(statenew(:,1),statenew(:,2),statenew(:,3))
xlabel('x [km]')
ylabel('y [km]')
ylabel('z [km]')

%% Phase 1: 100km - 40km
% Hop, then football [20km-40km]

r0 = [100 0 0];
v0 = [0 0 0];
x2  = 70;



y0_dot = n*x2/4;





%% Phase 2: 30km - 1km
% Hop


%% Phase 3: 1km - 300m
% Hop

%% Phase 4: 300m - 20m
% Hop


%% Phase 5: 20m - Approach
% Hop within 20m [+-1], vbar approach to target