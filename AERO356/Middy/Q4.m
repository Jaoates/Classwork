%% Question 4
clear all
close all
clc
%% Q4
T = 3800;%K
r = 6.26e8; %m
R = 6.704e10;%m

h = 6.62607015e-34;
k = 1.380649e-23;
c = 299792458;

% sun vals
% r= 695700e3
% R = 1.495978707e11
% T = 5800 


sig = 5.6704e-8;

A = 4*pi*r^2;

% exitance
syms dA
M = sig*T^4;
Q_dot = double(int(M,dA,[0,A]));
disp("The exitance of the star is "+string(Q_dot)+" W")
% irradiance
E_planet = M*(r^2/R^2);
disp("The irradiance at the planet is "+ string(E_planet)+" W/m^2")

% wavelength max
lam_peak = 2898/T; % micrometers

disp("The primary wavelength present as calculated from the formula is "+string(lam_peak)+" micrometers")
disp("By visual inspection of black body charts, the primary wavelength present is .8 micrometers")

% what star color
disp("This wavelength corresponds to deep red color. The star then would be a red star.")


% plancks law

h = 6.62607015e-34;
k = 1.380649e-23;
c = 299792458;

lam = linspace(1e-7,50e-7);

c1 = 2*pi*h*c^2;
c2 = h*c/k;
M = @(lam) (c1./lam.^5)./(exp(c2./(lam.*T))-1);
B = @(lam) (2*h*c^2./lam.^5)./(exp(h*c./(lam*k*T))-1);

figure
plot(lam,B(lam))
set(gca,"XScale","log");
set(gca,"YScale","log");
xlabel("Wave Length (m)")
ylabel("Spectral Radiance (W/(m^2*sr*Hz))")
title("Spectral Radiance vs Wavelength")

figure
plot(lam,M(lam))
xlabel("Wave Length (m)")
ylabel("Spectral Radiance (W/(m^2*Hz))")
title("Spectral Radiant Exitance vs Wavelength")




