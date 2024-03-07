clear all
close all
clc

% This is the caller for torque free motion in simulink

% add resources
% addpath("C:\joshFunctionsMatlab\")
addpath("..\")
load("MOInominal.mat")
load("initialOrientation.mat")
load("params.mat")

% initialize k
kd = params.kd;
kp = params.kp;
etac = params.etac;
epsc = params.epsc;

I = diag(nominal.moi);

% initialize initial conditions
w0 = [.001;-.001;.002];
E0 = [initial.t1i0;initial.t2i0;initial.t3i0];
C0 = initial.Ci0;
[eta0,eps0] = joshRotM2Quat(C0);

r0 = initial.r0;
v0 = initial.v0;


X0 = [w0;E0;eps0;eta0;r0;v0];
tspan = [0,30000];
% tspan = tspan(1):100:tspan(2);
options = odeset('RelTol', 1e-6,'AbsTol',1e-6);

ODEW = @(t,X) SC_ODE(t,X,I,params,nominal);
[t,X] = ode45(ODEW,tspan,X0,options);

% X = X(421:end,:) %%%%%%%% hack off first part for more clarity for now
% t = t(421:end)


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
title("Euler angles vs time")
xlabel("t(s)")
ylabel("degrees")
legend("x","y","z")

% plot Euler angles and convert from quat2angle
figure
hold on
q = [X(:,10),X(:,7),X(:,8),X(:,9)];
[x,y,z] = quat2angle(q,"XYZ");

plot(t,rad2deg(x),t,rad2deg(y),t,rad2deg(z))
title("Euler angles vs time q2a")
xlabel("t(s)")
ylabel("degrees")
legend("x","y","z")

% plot quaternion
figure
hold on 
plot(t,X(:,7),t,X(:,8),t,X(:,9),t,X(:,10))
title("Quaternion vs time")
xlabel("t(s)")
legend("i","j","k","q")

% % plot torque
% figure
% hold on 
% plot(t,X(:,17),t,X(:,18),t,X(:,19))
% title("Torque vs time")
% xlabel("t(s)")
% ylabel("Nm")
% legend("x","y","z")

% plot r
figure
hold on 
plot3(X(:,11),X(:,12),X(:,13))
title("r trajectory")
xlabel("x [km]")
ylabel("y [km]")
zlabel("z [km]")

[n,m] = size(X);
g = 6;
X = [X,zeros(n,g)];
for i = 1:n

    T = gravgrad(X(i,:)',I);
    X(i,m+1) = T(1);
    X(i,m+2) = T(2);
    X(i,m+3) = T(3);

    T = MagTorque(t(i),X(i,:)',nominal);
    X(i,m+4) = T(1);
    X(i,m+5) = T(2);
    X(i,m+6) = T(3);

    T = SolarPressureTorque(X(i,:)',nominal);
    X(i,m+7) = T(1);
    X(i,m+8) = T(2);
    X(i,m+9) = T(3);

    T = AeroTorque(X(i,:)',nominal);
    X(i,m+10) = T(1);
    X(i,m+11) = T(2);
    X(i,m+12) = T(3);

end

figure
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
xlabel("time(s)")
ylabel("torque (N*m)")
title("Torque From Gravity")
legend("x","y","z")

m = m+3;
figure
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
xlabel("time(s)")
ylabel("torque (N*m)")
title("Torque From Magnetic Dipole")
legend("x","y","z")

m = m+3;
figure
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
xlabel("time(s)")
ylabel("torque (N*m)")
title("Torque From SRP")
legend("x","y","z")

m = m+3;
figure
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
xlabel("time(s)")
ylabel("torque (N*m)")
title("Torque From Aero Drag")
legend("x","y","z")




