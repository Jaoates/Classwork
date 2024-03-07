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
% etac = params.etac;
% epsc = params.epsc;

I = diag(nominal.moi);

% initialize initial conditions
w0 = [.001;-.001;.002];
E0 = [initial.t1i0;initial.t2i0;initial.t3i0];
C0 = initial.Ci0;
[eta0,eps0] = joshRotM2Quat(C0);
eps0 = -eps0;
r0 = initial.r0;
v0 = initial.v0;

Ohm0 = [0;0;0];

X0 = [w0;E0;eps0;eta0;r0;v0;Ohm0];
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
title("Omega vs time b-ECI")
xlabel("time [s]")
ylabel("rads/s")
legend("x","y","z")
ylim([-2 2]*1e-3)

% plot Euler angles and convert from rads to degs
figure
hold on 
plot(t,rad2deg(X(:,4)),t,rad2deg(X(:,5)),t,rad2deg(X(:,6)))
title("Euler angles vs time b-ECI")
xlabel("time [s]")
ylabel("degrees")
legend("x","y","z")

% plot Euler angles and convert from quat2angle
figure
hold on
q = [X(:,10),X(:,7),X(:,8),X(:,9)];
[x,y,z] = quat2angle(q,"XYZ");

plot(t,rad2deg(x),t,rad2deg(y),t,rad2deg(z))
title("Euler angles vs time q2a b-ECI")
xlabel("time [s]")
ylabel("degrees")
legend("x","y","z")

% plot quaternion
figure
hold on 
plot(t,X(:,7),t,X(:,8),t,X(:,9),t,X(:,10))
title("Quaternion vs time b-ECI")
xlabel("time [s]")
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
title("r trajectory ECI")
xlabel("x [km]")
ylabel("y [km]")
zlabel("z [km]")

% plot Ohm
figure
hold on 
plot(t,X(:,17),t,X(:,18),t,X(:,19))
title("Reaction Wheel Speed")
legend("x","y","z")
ylabel("speed [rad/s]")
xlabel("time [s]")
ylim([-1.5 1.5])

[n,m] = size(X);
g = 1;
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


    r = X(i,11:13)';
    v = X(i,14:16)';
    qLE = RV2LVLH(r,v);
    qbE = [X(i,10),X(i,7:9)];
    qe = quatmultiply(quatconj(qLE), qbE);

    X(i,m+13) = qe(2);
    X(i,m+14) = qe(3);
    X(i,m+15) = qe(4);
    X(i,m+16) = qe(1);
    
    wLE = cross(r,v)/norm(r)^2;
    w = X(1:3);
    wbE = w;
    wbL = wbE-quatrotate(qbE,wLE')';
    X(i,m+17) = wbL(1);
    X(i,m+18) = wbL(2);
    X(i,m+19) = wbL(3);

    Mc = FEEDBACK_ODE(X(i,:)',params);
    X(i,m+20) = Mc(1);
    X(i,m+21) = Mc(2);
    X(i,m+22) = Mc(3);

end

figure
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
xlabel("time(s)")
ylabel("torque (N*m)")
title("Torque From Gravity b")
legend("x","y","z")
ylim([-1.5 1.5]*1e-8)

m = m+3;
figure
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
xlabel("time(s)")
ylabel("torque (N*m)")
title("Torque From Magnetic Dipole b")
legend("x","y","z")

m = m+3;
figure
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
xlabel("time [s]")
ylabel("torque (N*m)")
title("Torque From SRP b")
legend("x","y","z")

m = m+3;
figure
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
xlabel("time [s]")
ylabel("torque (N*m)")
title("Torque From Aero Drag b")
legend("x","y","z")

m = m+3;
figure
hold on 
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3),t,X(:,m+4))
title("Quaternion vs time b-LVLH")
xlabel("time [s]")
legend("i","j","k","q")

% plot Euler angles and convert from quat2angle
figure
hold on
q = [X(:,m+2),X(:,m+3),X(:,m+4),X(:,m+1)];
[x,y,z] = quat2angle(q,"XYZ");

plot(t,rad2deg(x),t,rad2deg(y),t,rad2deg(z))
title("Euler angles vs time q2a b-LVLH")
xlabel("time [s]")
ylabel("degrees")
legend("x","y","z")

m = m+4;
% plot omega
figure 
hold on
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
title("Omega vs time b-LVLH")
xlabel("time [s]")
ylabel("rads/s")
legend("x","y","z")
ylim([-2 2]*1e-3)

m = m+3;
figure 
hold on
plot(t,X(:,m+1),t,X(:,m+2),t,X(:,m+3))
title("Mc")
xlabel("time [s]")
ylabel("torque [N*m]")
legend("x","y","z")
ylim([-5 5]*1e-3)


quat = [t,X(:,7),X(:,8),X(:,9),X(:,10)];
save("quat.mat","quat","-mat")







