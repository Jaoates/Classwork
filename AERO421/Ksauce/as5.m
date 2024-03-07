clear; close all; clc

% addpath("../A320")


%% 1- Performance analysis
I = [1200 0 0; 0 2000 0; 0 0 2800];

zeta = .65;
t_s = 30; %s

%omega_n = log(t_s - .02*sqrt(1-zeta^2))/(-zeta*t_s);
omega_n = log(.02*sqrt(1-zeta^2))/(-zeta*t_s);

k_d = I * eye(3) * 2 * zeta * omega_n;
k_p = I * eye(3) * 2 * omega_n^2;

%% Model inputs

eps_0 = [.2; -.5; .3];
eta_0 = sqrt(1- eps_0'*eps_0);
omega_0 = [.1; -.05; .05];

e_0 = [0;0;0];


%% Case 1 
eps_c = [0;0;0];
eta_c = 1;
out = sim("FSFB.slx");
out_lin = sim("linearFB.slx");


%% Plots - Case 1
omega_labels = {'$\omega_x$', '$\omega_y$','$\omega_z$'};
% eul_labels = {'$\phi$', '$\theta$', '$\psi$'}; 
quat_labels = {'$q_1$', '$q_2$', '$q_3$', '$q_4$'};
T_labels = {'$T_x$', '$T_y$','$T_z$', 'norm(T)'};

omega = squeeze(out.omegaData);
q = squeeze(out.q);
T_c = squeeze(out.T_c);

omega_lin = squeeze(out_lin.omega);
eps_lin = squeeze(out_lin.eps);
tiledlayout(2,1)

nexttile
hold on
plot(out.tout, q)
plot(out_lin.tout, eps_lin, '--')
legend(quat_labels,  Interpreter="latex")
title("Nonlinear model (solid) vs linearized model (dash) ")

nexttile
hold on
plot(out.tout, omega)
plot(out_lin.tout, omega_lin, '--')
legend(omega_labels,  Interpreter="latex")
ylabel('$ \omega$ (rad/s)', Interpreter='latex')
xlabel('$t$ (s)', Interpreter='latex')

figure
hold on
plot(out.tout, T_c)
plot(out.tout, vecnorm(T_c))
legend(T_labels, Interpreter="latex")
title("Thruster Torque")
ylabel('$T_c$ (Nm)', Interpreter='latex')

disp('Both the nonlinear and linearized systems meet performance requirements for this case, settling within the required 30 seconds.')
%% Case 2
eps_c = [-.2; -.4; .2];
eta_c = sqrt(1- eps_c'*eps_c);
out = sim("FSFB.slx");
out_lin = sim("linearFB.slx");

%% Plots - Case 2

omega = squeeze(out.omegaData);
q = squeeze(out.q);
T_c = squeeze(out.T_c);

omega_lin = squeeze(out_lin.omega);
eps_lin = squeeze(out_lin.eps);

figure
tiledlayout(2,1)

nexttile
hold on
plot(out.tout, q)
plot(out_lin.tout, eps_lin, '--')
legend(quat_labels,  Interpreter="latex")
title("Nonlinear model (solid) vs linearized model (dash) ")

nexttile
hold on
plot(out.tout, omega)
plot(out_lin.tout, omega_lin, '--')
legend(omega_labels,  Interpreter="latex")
ylabel('$ \omega$ (rad/s)', Interpreter='latex')
xlabel('$t$ (s)', Interpreter='latex')


figure
hold on
plot(out.tout, T_c)
plot(out.tout, vecnorm(T_c))
legend(T_labels, Interpreter="latex")
title("Thruster Torque")
ylabel('$T_c$ (Nm)', Interpreter='latex')
xlabel('$t$ (s)', Interpreter='latex')

disp(['In this case, the nonlinear model does not meet the performance requirements. By 60 seconds, ' ...
    'it is still oscillating past 2% of the desired states for both the quaternion and angular velocity. ' ...
    'Theis is likely because the selected gains were derived from the linearized model, which assumed small ' ...
    'deviations from the equilibrium point of x = 0. As for the linearized model itself, it settled at the ' ...
    'appropriate time but due to the approximations made in calcuating the quaternions it settled only ' ...
    'within about 10% of the desired states. '])

%% 3- Thruster Selection

maxT = max(vecnorm(T_c)); %Nm
d = 2; %m
FS = 1.2;
F_req = FS * maxT/d; %N

fprintf("Max required force ~ %.2f N\n", F_req)
disp("One commercially available option that meets this requirement is the Moog " + ...
    "MONARC-90HT hydrazine monopropellant thruster, which provides a steady state thrust " + ...
    "of 116 N. More information can be found at: ")
disp("https://www.moog.com/content/dam/moog/literature/sdg/space/propulsion/moog-MonopropellantThrusters-Datasheet.pdf")