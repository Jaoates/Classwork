clear
% close all
clc
% load("MassProperties.mat");
% addpath("../../A320")

load("MOInominal.mat")
mprops_NO.I = diag(nominal.moi);
mprops_NO.COM = [
    0
    0
    0.2344];

%% Spacecraft Area properties

a = 2; % cube side length
l_s = 1; % sensor length
a_s = .25; % sensor cube side length

l_p = 3; % solar panel length
w_p = 2; % width
t_p = .05; % thickness

A = [a^2 + 2*t_p*l_p + a_s*l_s, a^2 + a_s*l_s, a^2 + 2*l_p*w_p];
% projected area in body x, y, and z
A = [A A]; % same area in opposite directions

Ax = [a^2, 2*t_p*l_p, a_s*l_s];
z_b = [0, 0, 1.5];
rhoz = sum(Ax.*z_b)/sum(Ax);

rho = [[0; 0; rhoz]-mprops_NO.COM, [0; 0; rhoz]-mprops_NO.COM, [0; 0; 0]];
rho = [rho [[0; 0; rhoz]-mprops_NO.COM, [0; 0; rhoz]-mprops_NO.COM, [0; 0; 0]]];

%% Model inputs

I = mprops_NO.I;
omega_0 = [.001; -.001; .002];
eps_0b = [0; 0; 0]; % body LVLH
eta_0b = 1; % body LVLH

Re = 6378; %km
mu_e = 398600;
orbit.h = 53335.2; %km^2/s
orbit.e = 0;
orbit.Omega = 0;
orbit.i = 98.43;
orbit.AOP = 0;
orbit.theta = 0;
orbit.a = orbit.h^2/mu_e * (1-orbit.e^2);

jd_0 = juliandate(2023, 1, 1, 12, 0, 0); % equinox?
s_eci = [1; 0; 0];

p = 4.5e-6; %N/m^2
[~, rho_a] = atmosnrlmsise00((orbit.a-Re)*1e3, 0, 0, 2023, 1, 1, 0);
rho_a = rho_a(6);

m = [0; 0; .5];

C_LVLH_ECI_0 = [0   -0.1466    0.9892; 0    0.9892    0.1466; -1.0000         0         0];

%% rotation sequence from body to ECI
[qtemp] = dcm2quat(C_LVLH_ECI_0);

% [eps_E0, eta_E0] = 
eta_E0 = qtemp(1) 
eps_E0 = qtemp(2:4)'


eta_Eb = eta_0b * eta_E0 - eps_E0' * eps_0b;
eps_Eb = eta_E0 * eps_0b + eta_0b * eps_E0 + joshCross(eps_0b)*eps_E0;

[phi, theta, psi] = quat2angle([eta_Eb eps_Eb'], 'ZYX');

e_0 = [psi; theta; phi]; %initial euler angles


%% Run simulation
out = sim("DisturbanceModels.slx");


%% plots

omega_labels = {'$\omega_x$', '$\omega_y$','$\omega_z$'};
eul_labels = {'$\phi$', '$\theta$', '$\psi$'};
quat_labels = {'$q_1$', '$q_2$', '$q_3$', '$q_4$'};
T_labels = {'$T_x$', '$T_y$','$T_z$'};

x = squeeze(out.x);
omega = x(5:end, :);
q = x(1:4, :);
% eul = squeeze(out.eul);
T_gg = squeeze(out.T_gg);
T_s = squeeze(out.T_s);
T_a = squeeze(out.T_a);

tiledlayout(3,1)
nexttile

plot(out.tout, omega)
legend(omega_labels,  Interpreter="latex")


nexttile
hold on
plot(out.tout, q)
legend(quat_labels,  Interpreter="latex")

% nexttile
% plot(out.tout, rad2deg(eul))
% legend(eul_labels, Interpreter="latex")

figure
plot(out.tout, T_gg)

figure
plot(out.tout, T_s)

figure
plot(out.tout, T_a)