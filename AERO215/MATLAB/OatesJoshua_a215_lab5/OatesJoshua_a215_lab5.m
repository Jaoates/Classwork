%% Joshua Oates - a215 - lab 5 - File Operations
clear all;
close all;
clc;

%% Part 1
% In a script, read into MATLAB the data from turbine experiment using 
% readtable.
TurbineDatIn = readtable('a215_lab5files\a215_lab5files\turbineRunData.xlsx');


% Use the . notation to pull the variables out of the table as 
% vectors.
dTime = TurbineDatIn.DTime_seconds_; % seconds

T01 = TurbineDatIn.T01_F_; % Fareinheit
T03 = TurbineDatIn.T03_F_;
T04 = TurbineDatIn.T04_F_;
T05 = TurbineDatIn.T05_F_;

P1 = TurbineDatIn.P1_psig_; % psi
P01 = TurbineDatIn.P01_psig_;
P03 = TurbineDatIn.P03_psig_;
P07 = TurbineDatIn.P07_psig_;

%Convert the temperatures T1, T03, T04, and T05 from F to K and the
% pressures from PSI to N/m2 using vector operations. 
% (temp in °F − 32) × 5/9 + 273.15 = temp in K
PSIconv = 6894.76; % 1 PSI = 6894.76 N/m2

P1 = P1 * PSIconv; % overwrite P in psi to P in N/m^2
P01 = P01 * PSIconv;
P03 = P03 * PSIconv;
P07 = P07 * PSIconv;

T01 = ( T01 - 32 ) * (5/9) + 273.15; % overwrite T in F to T in K
T03 = ( T03 - 32 ) * (5/9) + 273.15;
T04 = ( T04 - 32 ) * (5/9) + 273.15;
T05 = ( T05 - 32 ) * (5/9) + 273.15;

% Create two figures, one of all the temperatures in K and one of all the
% pressures in N/m2  versus "D Time"
figure (1)
hold ('on') 
plot (dTime, T01, dTime, T03, dTime, T04, dTime, T05);
xlabel ("delta time (s)");
ylabel ("temp (K)");
title ("Gas Turbine Temperatures");
legend("T01","T03","T04","T05");

figure (2)
hold ('on') 
plot (dTime, P1, dTime, P01, dTime, P03, dTime, P07);
xlabel ("delta time (s)");
ylabel ("pressure N/m^3");
title ("Gas Turbine Pressure");
legend("P1","P01","P03","P07");

%% Part 2
% Using Publish, submit this Lab 5 as a PDF
% lastnameFirstname_a215_lab5.pdf