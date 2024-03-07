%% HW5 joshua oates

% %% init parameters
% clear all
% kd = -.5;
% params.kd = kd;
% save("params.mat","params","-mat")


%% run
clear all;
% close all;
clc

% addpath("C:\joshFunctionsMatlab\")
addpath("C:\AERO421\FinalProj\")
load("MOIdetumble.mat")
load("initialOrientation.mat")
load("params.mat")

% initialize vars SC_ODE needs
I = diag(detumbleMOI.moi);
k = params.kd;



w0 = [-.05;.03;.2];
E0 = [initial.t1i0;initial.t2i0;initial.t3i0];
C0 = initial.Ci0;
[eta0,eps0] = joshRotM2Quat(C0);

tspan=[0,6000];

X0 = [w0;E0;eps0;eta0];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);

ODEW = @(t,X) SC_ODE(t,X,I,k);
[t,X] = ode45(ODEW,tspan,X0,options);

%%%%%%%%%%%%%%%%%%%%

figure 
hold on
plot(t,X(:,1),t,X(:,2),t,X(:,3))
title("Omega vs time (ODE45)")
xlabel("t")
ylabel("rads/s")
legend("x","y","z")

figure
hold on 
plot(t,rad2deg(X(:,4)),t,rad2deg(X(:,5)),t,rad2deg(X(:,6)))
title("Euler angles vs time dir (ODE45)")
xlabel("t")
ylabel("degrees")
legend("x","y","z")

% plot Euler angles and convert from quat2angle
figure
hold on
q = [X(:,10),X(:,7),X(:,8),X(:,9)];
[x,y,z] = quat2angle(q,"XYZ");

plot(t,x,t,y,t,z)
title("Euler angles vs time q2a (ODE45)")
xlabel("t(s)")
ylabel("degrees")
legend("x","y","z")

figure
hold on 
plot(t,X(:,7),t,X(:,8),t,X(:,9),t,X(:,10))
title("Quaternion vs time (ODE45)")
xlabel("t")
legend("i","j","k","q")








