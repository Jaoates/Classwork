% Joshua Oates - prelab 7 - ode solving

%% Section 0 - clean up
clear all;
close all;
clc

%% section 1 - create and plot solution to ode

% solve y" + (y)y' + y = 0
syms y(t) % use symbolic toolbox because it's alot easier
diffEq = diff(y,t,2) + y * diff(y,t) + y == 0; % create symbolic ODE
ode = matlabFunction(odeToVectorField(diffEq)); % create function handles that can be taken by normal matlab as normal ODE
initial = [0,1]; % use given intital conditions
timeSpan = [0,8]; % set domain to plot over
[T,Y] = ode45(@(T,Y) ode(Y), timeSpan, initial); % get output from ode45
figure
hold on
plot(T, Y(:,1)) % create and label plots
plot(T, Y(:,2))
legend("y","y'")
xlabel("t axis")
ylabel("y axis")
title("y"+'"'+" + y*y' + y = 0")