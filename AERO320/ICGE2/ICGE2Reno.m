% close all; clear; clc
addpath("C:\joshFunctionsMatlab\")
%% 1 
m_r = 600; %kg 
m_nc = 50;
m_w = 10; 
R = .25; 
R_w = .2;
t_w = .04;
h_r = 4;
h_nc = 1;

r_c_rocket = (2 * m_r + 4.25 * m_nc)/(m_r + m_nc);

I_r = (m_r) * diag([1/12*h_r^2 + .25*R^2, 1/12*h_r^2 + .25*R^2, .5*R^2]);

I_nc = m_nc * diag([3/80 * (h_nc^2 + 4*R^2), 3/80 * (h_nc^2 + 4*R^2), 3/10 * R^2]);

J_r = I_r - m_r*joshCross([0 0 2 - r_c_rocket]) *joshCross([0 0 2 - r_c_rocket]);
J_nc = I_nc - m_nc*joshCross([0 0 4.25 - r_c_rocket])*joshCross([0 0 4.25 - r_c_rocket]); 

% (a)
I_rocket = J_r + J_nc

I_w = (m_w) * diag([1/12*t_w^2 + .25*R_w^2, 1/12*t_w^2 + .25*R_w^2, .5*R_w^2])

%% 3

omegadotrel = [0; 0; .05];
omegarel = [0; 0; 0];
x0 = zeros(6,1);

R = @(phi, theta, psi) [cos(theta) sin(phi)*sin(theta) cos(phi)*sin(theta);...
                           0 cos(phi)*cos(theta) -sin(phi)*cos(theta); ... 
                           0 sin(phi) cos(phi)];

euldot = @(omega, phi, theta, psi) 1/cos(theta) * R(phi, theta, psi) * omega;



xdot = @(t,x, Td, omegadotrel, omegarel) [-(I_rocket+I_w)^-1 * (joshCross(x(1:3))* I_rocket * x(1:3) + I_w * omegadotrel + joshCross(x(1:3) + omegarel) * I_w  * (x(1:3) + omegarel) + Td); euldot(x(1:3), x(4), x(5), x(6))];
tspan = [0, 100];

[t,x] = ode45(@(t,x) xdot(t,x,0, omegadotrel, omegarel), tspan, x0);

figure
tl = tiledlayout('flow');
title(tl,"Problem 3")

nexttile
plot(t, x(:, 1), t, x(:,2), t, x(:,3))
xlabel("Time (s)")
ylabel("\omega (rad/s)")
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$'}, Interpreter="latex")

nexttile
plot(t, x(:, 4), t, x(:,5), t, x(:,6))
xlabel("Time (s)")
ylabel("\theta (rad)")
legend({'$\theta_x$', '$\theta_y$', '$\theta_z$'}, Interpreter="latex")

%% 4 

Td = [.1; 0; 0];
tspan = [0, 1000];

% (a)
x0 = zeros(6,1);
omegadotrel = zeros(3,1);
omegarel = zeros(3,1);
[t,x] = ode45(@(t,x) xdot(t,x,Td, omegadotrel,omegarel), tspan, x0);

figure
tl = tiledlayout('flow');
title(tl,"Problem 4 (a)")

nexttile
plot(t, x(:, 1), t, x(:,2), t, x(:,3))
xlabel("Time (s)")
ylabel("\omega (rad/s)")
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$'}, Interpreter="latex")

nexttile
plot(t, x(:, 4), t, x(:,5), t, x(:,6))
xlabel("Time (s)")
ylabel("\theta (rad)")
legend({'$\theta_x$', '$\theta_y$', '$\theta_z$'}, Interpreter="latex")

% (b)
x0 = [0; 0; .1; 0; 0; 0];
[t,x] = ode45(@(t,x) xdot(t,x,Td, omegadotrel,omegarel), tspan, x0);

figure
tl = tiledlayout('flow');
title(tl,"Problem 4 (b)")

nexttile
plot(t, x(:, 1), t, x(:,2), t, x(:,3))
xlabel("Time (s)")
ylabel("\omega (rad/s)")
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$'}, Interpreter="latex")

nexttile
plot(t, x(:, 4), t, x(:,5), t, x(:,6))
xlabel("Time (s)")
ylabel("\theta (rad)")
legend({'$\theta_x$', '$\theta_y$', '$\theta_z$'}, Interpreter="latex")

% (c)
x0 = zeros(6,1);
omegarel = [0; 0; 100];

[t,x] = ode45(@(t,x) xdot(t,x,Td, omegadotrel,omegarel), tspan, x0);

figure
tl = tiledlayout('flow');
title(tl,"Problem 4 (c)")

nexttile
plot(t, x(:, 1), t, x(:,2), t, x(:,3))
xlabel("Time (s)")
ylabel("\omega (rad/s)")
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$'}, Interpreter="latex")

nexttile
plot(t, x(:, 4), t, x(:,5), t, x(:,6))
xlabel("Time (s)")
ylabel("\theta (rad)")
legend({'$\theta_x$', '$\theta_y$', '$\theta_z$'}, Interpreter="latex")

% (d)
x0 = [0; 0; .1; 0; 0; 0];
[t,x] = ode45(@(t,x) xdot(t,x,Td, omegadotrel,omegarel), tspan, x0);

figure
tl = tiledlayout('flow');
title(tl,"Problem 4 (d)")

nexttile
plot(t, x(:, 1), t, x(:,2), t, x(:,3))
xlabel("Time (s)")
ylabel("\omega (rad/s)")
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$'}, Interpreter="latex")

nexttile
plot(t, x(:, 4), t, x(:,5), t, x(:,6))
xlabel("Time (s)")
ylabel("\theta (rad)")
legend({'$\theta_x$', '$\theta_y$', '$\theta_z$'}, Interpreter="latex")