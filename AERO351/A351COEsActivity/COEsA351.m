%% Joshua Oates - a_215 - HW 3 - COEs 

clear all;
close all;
clc;

addpath('C:\joshFunctionsMatlab')

R = [9031.5 -5316.9 -1647.2];
V = [-2.8640 5.1112 -5.0805];

% Write a script to call your function with the 𝑅⃗ and 𝑉⃗ vectors 
% from above as inputs. The script should print all six COE’s (with units) 
% neatly to the Command Window.

%create a template for a structure that can hold the COEs for a given orbit
COEs.a = 0;
COEs.e = 0;
COEs.nu = 0;
COEs.i = 0;
COEs.raan = 0;
COEs.aop = 0;
COEs.T = 0;

mu_e = 398600;

COEs1 = COEs; % use template to create COEs1. COEs1 is for the initial orbit
[COEs1.a,COEs1.e,COEs1.nu,COEs1.i,COEs1.raan,COEs1.aop,COEs1.hm,COEs1.T] = COEsOatesJoshua(R,V,mu_e);

disp("Semi-major axis       : " + COEs1.a + " km")
disp("eccentricity          : " + COEs1.e + " km")
disp("true anomaly          : " + COEs1.nu + " degrees")
disp("inclination           : " + COEs1.i + " degrees")
disp("RAAN                  : " + COEs1.raan + " degrees")
disp("Argument of periapsis : " + COEs1.aop + " degrees")
disp("Angular momentum      : " + COEs1.hm + " km^2/s")