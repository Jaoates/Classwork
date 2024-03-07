%% Midterm 1 - Joshua Oates

%% Q4
T = 3800;%K
r = 6.26e8; %m
R = 6.704e10;%m

% sun vals
% r= 695700e3
% R = 1.495978707e11
% T = 5800 


sig = 5.6704e-8;

A = 4*pi*r^2;

% exitance
syms dA
M = sig*T^4;
Q_dot = double(int(Q_unitA,dA,[0,A]));

% irradiance
E_planet = M*(r^2/R^2);

% wavelength max

lam_peak = b/T;

% what star color





