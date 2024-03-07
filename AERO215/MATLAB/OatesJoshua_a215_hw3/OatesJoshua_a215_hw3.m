%% Joshua Oates - a_215 - HW 3 - COEs 

clear all;
close all;
clc;

% 𝑅⃗ = -2315.9 𝐼̂ + 2168.6 𝐽̂ + 6314.5 𝐾̂ [𝑘𝑚]
% 𝑉⃗ = -3.0599𝐼̂  6.0645 𝐽̂ -3.2044 𝐾̂ [𝑘𝑚/𝑠]

R = [-2315.9 , 2168.6 , 6314.5];
V = [-3.0599 , 6.0645 , -3.2044];



% Write a MATLAB function to calculate the classical orbital 
% elements from any state vector (R and V pair) of a spacecraft in an orbit
% around Earth. The function should have 𝑅⃗ and 𝑉⃗ as inputs (expressed in the
% ECI frame) and output semi-major axis (𝑎), eccentricity (𝑒), true anomaly
% (nu), inclination (𝑖), RAAN (Ω), and argument of perigee (𝜔).

%% Part 4
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

COEs1 = COEs; % use template to create COEs1. COEs1 is for the initial orbit
[COEs1.a,COEs1.e,COEs1.nu,COEs1.i,COEs1.raan,COEs1.aop,COEs1.T] = COEsOatesJoshua(R,V);

disp("Semi-major axis       : " + COEs1.a + " km")
disp("eccentricity          : " + COEs1.e + " km")
disp("true anomaly          : " + COEs1.nu + " degrees")
disp("inclination           : " + COEs1.i + " degrees")
disp("RAAN                  : " + COEs1.raan + " degrees")
disp("Argument of periapsis : " + COEs1.aop + " degrees")