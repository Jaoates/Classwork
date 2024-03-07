%% HW5 joshua oates
clear all;
close all;
clc

addpath("C:\joshFunctionsMatlab\")

%% Cylinder
clear all
syms z r t R h m

[Cx,Cy,Cz] = joshAxisRotation();
% Integral of J(1,1) on paper included

x = r*cos(t) % to convert cartesien to polar for integration
y = r*sin(t)

rho = [x;y;z] % cartesian rho

rhox = joshCross(rho)
rhoxx = simplify(rhox*rhox)
J1 = simplify(int(int(int(-rhoxx*r,r,[0,R]),t,[-pi,pi]),z,[-h,0])) % triple integral over area, inertia matrix Cylinder

J1 = subs(J1,R^2*pi*h,m);% where density is the constant 1

J2 = limit(J1,h,0) % as h approaches 0, disk
J3 = limit(J1,R,0) % as R approaches 0, rod
disp("------------P1------------")
disp("My work for this problem have the following results: ")
disp("J for cylinder: ")
disp(J1)
disp("J for disk: ")
disp(J2)
disp("J for thin rod: ")
disp(J3)
disp("Computation of J11 can be found in included hand work.")
disp("Due to previous comments from the graders, all symbolic steps will show their outputs.")
%% cone
clear all
syms r t z m h R0

rc = [0;0;-(3/4)*h] % from hand calcs

x = r*cos(t) % to convert cartesien to polar for integration
y = r*sin(t)

rho = [x;y;z] % cartesian rho
rhox = joshCross(rho)
rhoxx = simplify(rhox*rhox)

R = -R0*z/h
J = simplify(int(int(int(-rhoxx*r,r,[0,R]),t,[-pi,pi]),z,[-h,0])) % triple integral over area, inertia matrix cone

w = [1;t;sin(t)]
dw = diff(w)
wx = joshCross(w)

rcx = joshCross(rc)
I = J + m*rcx*rcx

Tc = I*dw + wx*I*w
Tc = subs(Tc,m,1)
Tc = subs(Tc,h,1)
Tc = subs(Tc,R0,1)
Tc = simplify(Tc)

disp("------------P2------------")
disp("My work for this problem have the following results: ")
disp("Hand calculations are included which find the center of mass of a cone along with an initial guess.")
disp("J cone: ")
disp(J)
disp("I cone: ")
disp(I)
disp("Net torque: ")
disp(Tc)

disp("Due to previous comments from the graders, all symbolic steps will show their outputs.")


%% ODE
clear all
close all

m = 1; 
h = 3;
r = 1;
w0 = [.5;-1;.5];
E0 = [0;0;0];
C0 = eye(3);
[eta0,eps0] = joshRotM2Quat(C0);


I = [[1 0 0];...
    [0 1 0];...
    [0 0 0]]*(1/12)*m*(3*r^2+h^2);
I(3,3) = .5*m*r^2;

tspan=[0,15];

X0 = [w0;E0;eps0;eta0];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
[tC,XC] = ode45(@odefunCoast,tspan,X0,options);
[tT,XT] = ode45(@odefunTorque,tspan,X0,options);

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
plot(tC,XC(:,4),tC,XC(:,5),tC,XC(:,6))
title("Euler angles vs time (no torque)")
xlabel("t")
ylabel("rads")
legend("x","y","z")

figure
hold on 
plot(tC,XC(:,7),tC,XC(:,8),tC,XC(:,9))
title("Epsilon vs time (no torque)")
xlabel("t")
legend("x","y","z")

figure
plot(tC,XC(:,10))
title("Eta vs time (no torque)")
xlabel("t")

%%%%%%%%%%%%%%%%%%%%%

figure 
hold on
plot(tT,XT(:,1),tT,XT(:,2),tT,XT(:,3))
title("Omega vs time (torque)")
xlabel("t")
ylabel("rads/s")
legend("x","y","z")

figure
hold on 
plot(tT,XT(:,4),tT,XT(:,5),tT,XT(:,6))
title("Euler angles vs time (torque)")
xlabel("t")
ylabel("rads")
legend("x","y","z")

figure
hold on 
plot(tT,XT(:,7),tT,XT(:,8),tT,XT(:,9))
title("Epsilon vs time (torque)")
xlabel("t")
legend("x","y","z")

figure
plot(tT,XT(:,10))
title("Eta vs time (torque)")
xlabel("t")

disp("------------P3------------")
disp("My work for this problem have the following results: ")
disp("See included hand calculations for equivalent cuboid.")
disp("See the 8 included plots.")


%% dependancies
depends = matlab.codetools.requiredFilesAndProducts("C:\AERO320\HW5\HW5.m");
depends = depends'


%% functions

function Xdot = odefunCoast(t,X)
w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);
T = [0;0;0];

I = [[1 0 0];[0,1,0];[0,0,.5]];

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

function Xdot = odefunTorque(t,X)
w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);
T = [-1;0;.5];

I = [[1 0 0];[0,1,0];[0,0,.5]];

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



