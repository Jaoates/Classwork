% Joshua Oates - Aero 300 - lab 7

%% section 0 - clean up
close all;
clear all;
clc;

%% section 1 - use rkf45 and ode45 on ode

% y(0)=1 
y0 = 1; % use given intital conditions
tspan = [0,pi/3]; % set domain to plot over
rTol = 1e-6;
opts = odeset('AbsTol',rTol);
[T1,Y1] = rkf45(@my_ode, tspan, y0,.1 ,rTol); % get output from rkf45
[T2,Y2] = ode45(@my_ode, tspan, y0, opts); % get output from ode45

figure
hold on
plot(T1, Y1) 
plot(T2, Y2) % create and label plots

legend("y-rkf45","y-ode45")
xlabel("t axis")
ylabel("y axis")
title("y' =(y-t-1)^2 +2")

y = @(t) 1+t+tan(t); % true solution
e1 = abs(Y1-y(T1));
e2 = abs(Y2-y(T2));

figure
hold on
plot(T1, e1) 
plot(T2, e2) % create and label plots
ax = gca();
ax.YScale = "log";

legend("y-rkf45","y-ode45")
xlabel("t axis")
ylabel("error")
title("error for: y' =(y-t-1)^2 +2")

disp("My function has aproximately 100 times more absolute error than ode45. Additionally, ode45 has a much less constant error than rfk45. This is likely becuase ODE 45 will increase h again after the relative error is sufficently low where as my rkf 45 will only ever decrease h.")

clear all

%% section 2 - Lorenz equation

y0 = [1,1,1];
tspan = [0,50];
rTol = 1e-4;

opts = odeset('AbsTol',1e-6);
[T1,Y1] = rkf45(@lorenz_ode,tspan,y0,.001,rTol);
[T2,Y2] = ode45(@lorenz_ode,tspan,y0,opts); % Runge-Kutta 4th/5th order ODE solver

figure
hold on
plot3(Y1(1,:),Y1(2,:),Y1(3,:))
plot3(Y2(:,1),Y2(:,2),Y2(:,3))
axis("equal")
campos([30,-25,40])
camtarget([1,0,0])

legend("y-rkf45","y-ode45")
xlabel("x axis")
ylabel("y axis")
zlabel("z axis")
title("Lorenz Equation")

%% - function def
function [dydt] = my_ode(t,Y)

y = Y(1);
yp = (y-t-1)^2 + 2;

dydt = yp;
end