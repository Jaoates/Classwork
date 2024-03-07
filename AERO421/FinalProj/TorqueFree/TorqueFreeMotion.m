%% HW5 joshua oates

% This file is the matlab version of doing the ODEs over time. Its great
% and quick, I wrote it in a few minutes

clear all;
close all;
clc

addpath("C:\joshFunctionsMatlab\")
addpath("C:\AERO421\FinalProj\")
load("MOI.mat")
load("initialOrientation.mat")

global ImatGlo
ImatGlo= diag(nominal.moi);


w0 = [.001;-.001;.002];
E0 = [initial.t1i0;initial.t2i0;initial.t3i0];
C0 = initial.Ci0;
[eta0,eps0] = joshRotM2Quat(C0);

tspan=[0,6000];

X0 = [w0;E0;eps0;eta0];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
[tC,XC] = ode45(@odefunCoast,tspan,X0,options);

%%%%%%%%%%%%%%%%%%%%

figure 
hold on
plot(tC,XC(:,1),tC,XC(:,2),tC,XC(:,3))
title("Omega vs time (no torque)")
xlabel("t")
ylabel("rads/s")
legend("x","y","z")

figure
hold on 
plot(tC,rad2deg(XC(:,4)),tC,rad2deg(XC(:,5)),tC,rad2deg(XC(:,6)))
title("Euler angles vs time (no torque)")
xlabel("t")
ylabel("degrees")
legend("x","y","z")

figure
hold on 
plot(tC,XC(:,7),tC,XC(:,8),tC,XC(:,9),tC,XC(:,10))
title("Quaternion vs time (no torque)")
xlabel("t")
legend("i","j","k","q")




%% functions

function Xdot = odefunCoast(t,X)
global ImatGlo
I = ImatGlo;

w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);
T = [0;0;0]; % torque

% this matrix is part of the EOM for a 321 sequence for Edot. It is not a
% rotation matrix
mat = [[cos(E(2)),      sin(E(2))*sin(E(1)),            sin(E(2))*cos(E(1))];...
    [0,                 cos(E(2))*cos(E(1)),           -cos(E(2))*sin(E(1))];...
    [0,                 sin(E(1)),                      cos(E(1))           ]]*(1/cos(E(2)));


Edot = mat*w;

epsx = joshCross(eps);
epsdot = .5*(eta*eye(3)+epsx)*w;
etadot = -.5*eps'*w;

wx = joshCross(w);
wdot = -inv(I)*(wx*I*w-T);

Xdot = [wdot;Edot;epsdot;etadot];
end





