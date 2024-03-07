%%Joshua Oates - Aero 215 - Lab 4

clear all;
close all;
clc;

%% Part 1
%write your own cross product code

H1 = [3,2,1];
H2 = [9,8,7];

JoshCross = crossProductOatesJoshua (H1,H2);
MatCross = cross (H1,H2);

%% Part 2
% A vector of measured velocity values V=[1,2,3,4,5] with units of ft/sec 
% corresponds to the time vector t=[0.2,0.4,0.6,0.9,1.1] with units of seconds.
% Plot V, V^2 , and V^3 versus time vector t on the same figure, each with a
% different color. Do not hard code (that is, don't do the math by hand and
% then write the values into MATLAB) the values for V^2 or V^3. Properly label 
% the axes, include a legend and title.

%Create vectors
V=[1,2,3,4,5]; % with units of ft/sec 
t=[0.2,0.4,0.6,0.9,1.1]; %with units of seconds

%element wise square and cube
V2 = V.^2; % ft^2/s^2
V3 = V.^3; % ft^3/s^3

%plot them
hold ("on");
plot (t,V);
plot (t,V2);
plot (t,V3);
xlabel("time (s)");
ylabel("V (ft/s), V^2 (ft^2/s^2), V^3 (ft^3/s^3)");
legend("V (ft/s)", "V^2 (ft^2/s^2)", "V^3 (ft^3/s^3)");  

 
 
