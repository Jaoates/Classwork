clear all
close all
clc

% This is the caller for torque free motion in simulink

% add resources
addpath("C:\joshFunctionsMatlab\")
addpath("C:\AERO421\FinalProj\")
load("MOI.mat")
load("initialOrientation.mat")

% initialize I
I = diag(nominal.moi);
T = [0;0;0];

% initialize initial conditions
w0 = [.001;-.001;.002];
E0 = [initial.t1i0;initial.t2i0;initial.t3i0];
C0 = initial.Ci0;
[eta0,eps0] = joshRotM2Quat(C0);

X0 = [w0;E0;eps0;eta0];

% sim call
simOut = sim('torqueFree.slx');

% sim output
t=simOut.ScopeData{1}.Values.Time;
X=simOut.ScopeData{1}.Values.Data;

% plot omega
figure 
hold on
plot(t,X(:,1),t,X(:,2),t,X(:,3))
title("Omega vs time (no torque)")
xlabel("t")
ylabel("rads/s")
legend("x","y","z")

% plot Euler angles and convert from rads to degs
figure
hold on 
plot(t,rad2deg(X(:,4)),t,rad2deg(X(:,5)),t,rad2deg(X(:,6)))
title("Euler angles vs time (no torque)")
xlabel("t")
ylabel("degrees")
legend("x","y","z")

% plot quaternion
figure
hold on 
plot(t,X(:,7),t,X(:,8),t,X(:,9),t,X(:,10))
title("Quaternion vs time (no torque)")
xlabel("t")
legend("i","j","k","q")

