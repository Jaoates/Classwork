clear all
close all
clc

% This is the caller for torque free motion in simulink

% add resources
% addpath("C:\joshFunctionsMatlab\")
addpath("..\")
load("MOIdetumble.mat")
load("initialOrientation.mat")
load("params.mat")

% initialize I
I = diag(detumbleMOI.moi);

% initialize k
k = params.kd;

% initialize initial conditions
w0 = [-.05;.03;.2];
E0 = [initial.t1i0;initial.t2i0;initial.t3i0];
C0 = initial.Ci0;
[eta0,eps0] = joshRotM2Quat(C0);

X0 = [w0;E0;eps0;eta0];

% sim call
simOut = sim('detumble.slx');

% sim output
t=simOut.ScopeData{1}.Values.Time;
X=simOut.ScopeData{1}.Values.Data;
X=squeeze(X)';
% plot omega
figure 
hold on
plot(t,X(:,1),t,X(:,2),t,X(:,3))
title("Omega vs time")
xlabel("t(s)")
ylabel("rads/s")
legend("x","y","z")

% plot Euler angles and convert from rads to degs
figure
hold on 
plot(t,rad2deg(X(:,4)),t,rad2deg(X(:,5)),t,rad2deg(X(:,6)))
title("Euler angles vs time (Sim)")
xlabel("t(s)")
ylabel("degrees")
legend("x","y","z")

% plot Euler angles and convert from quat2angle
figure
hold on
q = [X(:,10),X(:,7),X(:,8),X(:,9)];
[x,y,z] = quat2angle(q,"XYZ");

plot(t,rad2deg(x),t,rad2deg(y),t,rad2deg(z))
title("Euler angles vs time q2a (Sim)")
xlabel("t(s)")
ylabel("degrees")
legend("x","y","z")

% plot quaternion
figure
hold on 
plot(t,X(:,7),t,X(:,8),t,X(:,9),t,X(:,10))
title("Quaternion vs time (Sim)")
xlabel("t(s)")
legend("i","j","k","q")

% plot torque
figure
hold on 
plot(t,X(:,11),t,X(:,12),t,X(:,13))
title("Torque vs time (Sim)")
xlabel("t(s)")
ylabel("Nm")
legend("x","y","z")

