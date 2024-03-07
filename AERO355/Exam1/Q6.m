%% Exam 1 - A355 - question 6 - Joshua Oates
clear all
close all
clc

d = load("AODensity.txt");
alt = 320;%km
d = d(d(:,1) == alt,:);
AOmin = d(2); % atoms/m^3
AOmax = d(3);
AOave = (AOmax+AOmin)/2;

r_e = 6378;%km
mu_e = 3.98600e14;%m^3/s^2
a = r_e+alt; % km
a = a*1000; % m
V = (mu_e/a)^.5; % m/s

AOn = [AOmin,AOave,AOmax];
clear d AOmin AOmax AOave

flux = AOn*V; % atoms/m^2/s
E = 10.5e-24; %cm^3/atom
E = E/(1e6); %m^3/atom

t = 365.25*24*60*60; % year in seconds
dx_dt = E*flux; % errosion rate in m/s
dx_dt = dx_dt.*100; % cm/s
dx_dt = dx_dt*t; % cm/yr
t = linspace(0,1);

%%%%%%%%%%%%%%%%%%
figure 
hold on
plot(t,dx_dt(1).*t)
plot(t,dx_dt(2).*t)
plot(t,dx_dt(3).*t)
legend('Solar Min','Solar Ave','Solar Max','Location','best')
xlabel('time [years]')
ylabel('depth [cm]')
title('Errosion Depth Over Time')

disp("The errosion depth for the sample at solar max after 1 year is: "+string(dx_dt(3))+" cm.")
disp("The errosion depth for the sample at solar min after 1 year is: "+string(dx_dt(1))+" cm.")
disp("The errosion will reach a depth of .04 cm, thereby breaking the electrical connection in the solar cell, at solar max after: "+string((.04/dx_dt(3))*365.25)+" days.")
disp("The errosion will reach a depth of .04 cm, thereby breaking the electrical connection in the solar cell, at solar min after: "+string((.04/dx_dt(1))*365.25)+" days.")
disp("It is expected that there is less time to repair the spacecraft if it is at solar max because the concentration of AO will be higher at solar max than at solar min. It is also expected that errosion rate will be higher at solar max for the same reason. This can be seen in the plot of errosion depth over time. Becuase errosion happens linearly and we do not assume that the errosion rate ever changes, it makes sense that errosion depth over time is linear.")

