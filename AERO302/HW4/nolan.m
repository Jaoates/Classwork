%% AERO 302 - Assignment 4 - Nolan Kulp

clc; close all; clear

%% SP

% 5
t = linspace(0,100,100)';

% Sx = (-10*exp(-t).*(10*sin(10.*t) - cos(10.*t))) / (101);
Sx = 100.*cos(10.*t).*exp(-t);
Sy = 0*t + 10;

figure(1)
plot(Sx,Sy)
xlabel("x-position (cm)")
ylabel("y-position (cm)")
title("Streamline")

% 6
Px = -10*exp(-10)*(cos(10*t) + 10*(sin(10*t)));
Py = 0*t + 10;

figure(2)
tl = tiledlayout(1,2);
title(tl,"Pathline")

nexttile;
plot(Px,Py)
xlabel("x-position (cm)")

ylabel("y-position (cm)")
title("Pathline, x vs. y")

nexttile;
plot(t,Px)
xlabel("time (s)")
ylabel("x-position (cm)")
title("Pathline, t vs. x")

figure(3)
plot3(t,Px,Py)
xlabel("time (s)")
ylabel("x-position (cm)")
zlabel("y-position (cm)")
title("Pathline, t vs. x vs. y")



