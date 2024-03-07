%% Aero303 - Joshua Oates
clear all
close all
clc

addpath('C:\AERO303\compressibleFlowRelations\')


%% Problem 1
%{
A M1 = 2.8 flow is deflected by a compression corner with an angle of 15 deg. If the deflection angle 
is doubled to 30 deg, what is the change in shock strength? 
%}

M1 = 2.8;
theta1 = 15;
theta2 = 30;

oblique1 = joshComp([theta1 M1],6,'TMW');
Mn1 = M1*sind(oblique1.beta);
normal1 = joshComp(Mn1,2,'M');


oblique2 = joshComp([theta2 M1],6,'TMW');
Mn2 = M1*sind(oblique2.beta);
normal2 = joshComp(Mn2,2,'M');

disp("Problem 1")
disp("P0 loss for theta 15 is: "+string((1-normal1.P02_P01)*100)+' percent')
disp("P0 loss for theta 30 is: "+string((1-normal2.P02_P01)*100)+' percent')


%% Problem 2
%{
A M1 = 2.8 flow is deflected by a compression corner with an angle of 15 deg. If M1 is reduced to 2.0, 
what is the change in shock strength?  
%}

clear all
M1 = 2.8;
M2 = 2.0;
theta1 = 15;

oblique1 = joshComp([theta1 M1],6,'TMW');
Mn1 = M1*sind(oblique1.beta);
normal1 = joshComp(Mn1,2,'M');


oblique2 = joshComp([theta1 M2],6,'TMW');
Mn2 = M2*sind(oblique2.beta);
normal2 = joshComp(Mn2,2,'M');

disp("Problem 2")
disp("P0 loss for M = 2.8 is: "+string((1-normal1.P02_P01)*100)+' percent')
disp("P0 loss for M = 2.0 is: "+string((1-normal2.P02_P01)*100)+' percent')

%% Problem 3
%{
A supersonic 2D inlet is designed to operate at Mach 3.0. Compare the stagnation pressure 
downstream of the normal shock for the two design options: 
a. The compression and slowing down of the flow takes place through a single normal shock 
b. The compression and slowing down of the flow takes place across wedge-shaped ramps, that 
decelerates the flow through two weak oblique shocks followed by a normal shock. The 
ramps are 8 degrees each. 
%}
clear all
M1 = 3.0;
theta = 8;

normal = joshComp(M1,2,'M');

% first oblique in b
oblique1 = joshComp([theta M1],6,'TMW');
Mn1 = oblique1.M*sind(oblique1.beta);
normal1 = joshComp(Mn1,2,'M');

M2 = normal1.M2/sind(oblique1.beta-oblique1.theta);

% second oblique in b
oblique2 = joshComp([theta M2],6,'TMW');
Mn2 = oblique2.M*sind(oblique2.beta);

normal2 = joshComp(Mn2,2,'M');
M3 = normal2.M2/sind(oblique2.beta-oblique2.theta);

% additional normal in b
normal3 = joshComp(M3,2,'M');

ratioa = normal.P02_P01;
ratiob = normal1.P02_P01 * normal2.P02_P01 * normal3.P02_P01;

disp("Problem 3")
disp("P0 loss for 1 normal shock is: "+string((1-ratioa)*100)+" percent.")
disp("P0 loss for multiple shocks is: "+string((1-ratiob)*100)+" percent.")

%% Problem 4
%For a horizontal flow of air, M1 = 2.0, and incident shock angle βi = -40º, determine βr, M2 and M3. 

clear

M1 = 2.0;
beta1 = 40;

oblique1 = joshComp([beta1 M1],6,'BM');
Mn1 = oblique1.M*sind(oblique1.beta);
normal1 = joshComp(Mn1,2,'M');
M2 = normal1.M2/sind(oblique1.beta-oblique1.theta);

oblique2 = joshComp([oblique1.theta M2],6,'TMW');
Mn2 = oblique2.M*sind(oblique2.beta);
normal2 = joshComp(Mn2,2,'M');
M3 = normal2.M2/sind(oblique2.beta-oblique2.theta);

disp("Problem 4")
disp("Beta_r is: "+string(oblique2.beta))
disp("M2 is: "+string(M2))
disp("M3 is: "+string(M3))



