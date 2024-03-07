%% Exam 1 - A355 - question 7 - Joshua Oates
close all
clear all
clc


m = 10; % kg
A = .1;%m^2
TML = .25; %
TML = TML/100;
CVCM = .01; %
CVCM = CVCM/100;
t0 = 3e-12; %s

Ea = [5,10,15,20]'; % kcal/mol

T = linspace(50,300)';
R=1.9858775*10^-3; % kcal/K*mol



%%%%%%%%%%%%%%%%
figure
hold on

for i = 1:length(Ea)
    t(:,i) = t0*exp(Ea(i)./(R.*T));
    plot(T,t(:,i))
    xlabel('temp [K]')
    ylabel('residence time [s]')
    title('Residence Time Vs Temperature')
end
legend('Ea = ' + string(Ea) + ' kcal/mol')

gca().set('YScale','log')
gca().set('XScale','log')
disp("The residence time with an activation energy of 20 kcal/mol at 300 K is: "+string(t(end,end))+" s.")

%% part B

Ea = 15; % kcal/mol
T = [0,50,100]';
T = T+273.15;% K
m = 5; % g

t1 = 0:6;
t2 = 1:7;

dm = m*TML;

q0 = dm/(2*exp(-Ea/(R*(125+273.15))) *((24*3600)^.5));

x = t1.^.5-t2.^.5;

for j = 1:3
    for i = 1:7
        dm(i,j) = 2*q0*exp(-Ea/(R*T(j)))*x(i);
    end
end

figure
hold on
plot(t2,dm(:,1),'*')
plot(t2,dm(:,2),'*')
plot(t2,dm(:,3),'*')
gca().set('YScale','log')
legend('temp = '+string(T)+' K', 'location','best')
title('Daily Mass Loss for Samples at 3 Temps')
xlabel('Day')
ylabel('mass loss [g]')

disp("The mass loss on day 7 at 100 C is: "+string(dm(end,end))+" g.")
%% part c
clear dm
T = 250; % K
vf = linspace(.0003,.03,10);
Ea = 15;
rho = 1; % g/cm^3



for w = 1:length(vf)
    for i = 1:7
        dm = 2*q0*exp(-Ea/(R*T))*x(i);
        ra(i,w) = -vf(w)*dm*(1/rho);
    end
end

figure
hold on 

for w = 1:length(vf)
        plot(t2,ra(:,w),'*')
end

legend('view factor = '+string(vf))
xlabel('Day')
ylabel('arrival rate [cm^3]')
title('Daily Arrival Rate at Various View Factors at 250 K')
disp("The arrival rate on the 7th day at a view factor of .03 is: "+string(ra(end,end)) + "cm^3.")

%% discussion
disp("a. It is expected that residence time decreases dramatically with and increase temperature and will increase with greater activation energy. The plot for part a. indicates that this is the case.")
disp("b. It is expected that mass loss will decrease over time as more of the volitile material is off gassed. Day over day it is expected that mass loss will be lower. It is also expected that mass loss will be greater for higher temperatures as there is more energy to create off gassing. These trends can be seen in the figure for part b with daily mass loss at 3 temps. Note that all values are negative because I am measureing change in mass which is always going down, but the highest temperature sees much more mass loss than the other two and all mass loss at a given temperature is decreasing over time.")
disp("c. It is expected that daily arrival of material will be higher for higher view factors and decrease day over day. The arrival rate is seen to decrease bc arrival on a different surface comes from the the sample which is off gassing. As the sample off gasses less, the arrival rate will go down. I plotted ten view factors in the range that we are studying so that we can see the arrival rate increasing with better view factors.")