%% AERO 402 --- Homework 3 --- Nolan Kulp

clc; close all; clear

%% 2 - Nitrogen

% Ideal

k_nit = 1.4;
R_nit = 297; % J / kg*K
Tc_nit = 298; % K
F_nit = 40; % N
g_nit = 9.81;
I_tot_nit = 40000; % N*s
cstar_nit = (sqrt(k_nit*R_nit*Tc_nit)) / (k_nit * (sqrt((2/(k_nit+1))^((k_nit+1)/(k_nit-1)))));

p_ratio_nit = linspace(1,3000,1000)'; % pc/pe
p_ratio_flip_nit = 1 ./ p_ratio_nit; % pe/pc

CF_vector_nit = sqrt(   ((2*k_nit^2) / (k_nit-1)) * ( (2/(k_nit+1))^((k_nit+1)/(k_nit-1)) ) * (1 - (p_ratio_flip_nit.^((k_nit-1)/k_nit)))   );

% figure(1)
% plot(p_ratio_nit,CF_vector_nit)
% grid on
% xlabel("P_c/P_e")
% ylabel("C_F")

CF_ideal_nit = 1.8; % curve asymptote

mdot_ideal_nit = F_nit / (cstar_nit * CF_ideal_nit);

Isp_ideal_nit = F_nit / (mdot_ideal_nit * g_nit);

tb_nit = I_tot_nit / F_nit;

m_ideal_nit = tb_nit * mdot_ideal_nit;

T_amb = 298;

P_tank = linspace(5,200,1000)'; % atm

V_ideal_nit = (m_ideal_nit * R_nit * T_amb) ./ (P_tank/9.8692e-6);

% figure(3)
% plot(P_tank,V_ideal_nit)
% grid on
% xlabel("Tank Pressure [atm]")
% ylabel("Tank Volume [m^3]")

% Real

Eps_vector_nit = ( 2/(k_nit+1))^(1/(k_nit-1)) .* p_ratio_nit.^(1/k_nit) .* ( ((k_nit+1)./(k_nit-1)) .*(1-p_ratio_nit.^((1-k_nit)/k_nit)) ).^-.5;

% figure
% plot(Eps_vector_nit,CF_vector_nit)
% grid on
% xlabel("Area Ratio")
% ylabel("C_F")

Eps_real_nit = 25;

CF_real_nit = 1.67;

mdot_real_nit = F_nit / (cstar_nit * CF_real_nit);

mdot_vector_nit = F_nit ./ (cstar_nit * CF_vector_nit);

Isp_real_nit = F_nit / (mdot_real_nit * g_nit);

Isp_vector_nit = F_nit ./ (mdot_vector_nit * g_nit);

tb_nit = I_tot_nit / F_nit;

m_real_nit = tb_nit * mdot_real_nit;

V_real_nit = (m_real_nit * R_nit * T_amb) ./ (P_tank/9.8692e-6);

% Plot 3
figure
plot(P_tank,V_ideal_nit)
hold on
plot(P_tank,V_real_nit)
grid on
xlabel("Tank Pressure [atm]")
ylabel("Tank Volume [m^3]")
title("Nitrogen Tank Pressure vs. Volume")
legend("Ideal","Real")

At_nit = (mdot_real_nit .* cstar_nit) ./ (P_tank./9.8692e-6);
Ae_nit = Eps_real_nit * At_nit;

D_throat_nit = 2 .* (sqrt(At_nit./pi));
D_exit_nit = 2 .* (sqrt(Ae_nit./pi));

%% 2 - Hydrazine

% Ideal

k = 1.2;
R = 260; % J / kg*K
Tc = 1650; % K
F = 40; % N
g = 9.81;
rho = 1004.5; % kg/m3
I_tot = 40000; % N*s
cstar = (sqrt(k*R*Tc)) / (k * (sqrt((2/(k+1))^((k+1)/(k-1)))));

p_ratio = linspace(1,3000,1000)'; % pc/pe
p_ratio_flip = 1 ./ p_ratio; % pe/pc

CF_vector = sqrt(   ((2*k^2) / (k-1)) * ( (2/(k+1))^((k+1)/(k-1)) ) * (1 - (p_ratio_flip.^((k-1)/k)))   );

% Plot 1
figure
plot(p_ratio,CF_vector)
hold on
plot(p_ratio_nit,CF_vector_nit)
grid on
xlabel("P_c/P_e")
ylabel("C_F")
legend("Hydrazine","Nitrogen")
title("Pressure Ratio vs. CF")

CF_ideal = 2; % curve asymptote

mdot_ideal = F / (cstar * CF_ideal);

Isp_ideal = F / (mdot_ideal * g);

tb = I_tot / F;

m_ideal = tb * mdot_ideal;

T_amb = 298;

P_tank = linspace(5,200,1000)'; % atm

V_ideal = m_ideal / rho;

% Real

Eps_vector = ( 2/(k+1))^(1/(k-1)) .* p_ratio.^(1/k) .* ( ((k+1)./(k-1)) .*(1-p_ratio.^((1-k)/k)) ).^-.5;

Eps_real = 25;

CF_real = 1.7;

mdot_real = F / (cstar * CF_real);

mdot_vector = F ./ (cstar * CF_vector);

Isp_real = F / (mdot_real * g);

Isp_vector = F ./ (mdot_vector * g);

tb = I_tot / F;

m_real = tb * mdot_real;

% Plot 2
figure
plot(Eps_vector,Isp_vector)
hold on
plot(Eps_vector_nit,Isp_vector_nit)
grid on
xlabel("Expansion Ratio")
ylabel("I_s_p [s]")
legend("Hydrazine","Nitrogen")
title("Expansion Ratio vs. I_s_p")

V_real = m_real / rho;

At = (mdot_real .* cstar) ./ (P_tank./9.8692e-6);
Ae = Eps_real * At;

D_throat = 2 .* (sqrt(At./pi));
D_exit = 2 .* (sqrt(Ae./pi));

% Plot 4
figure
plot(P_tank,D_throat)
hold on
plot(P_tank,D_exit)
plot(P_tank,D_throat_nit)
plot(P_tank,D_exit_nit)
grid on
xlabel("Tank Pressure [atm]")
ylabel("Diameter [m]")
legend('Nitrogen (throat)','Nitrogen (exit)','Hydrazine (throat)','Hydrazine (exit)')
title("Tank Pressure vs. Nozzle Diameters")

% figure
% plot(Eps_vector,CF_vector)
% grid on
% xlabel("Area Ratio")
% ylabel("C_F")

%% 3

Eps = 20;
mdot = 0.05; % kg/s

% Hydrazine

gamma = 1.2;

k_hydrazine = 1.2;
cstar_hydrazine = cstar;

syms P_ratio_hydrazine

P_ce_ratio_hydrazine = vpasolve(Eps == ( 2/(gamma+1))^(1/(gamma-1)) * P_ratio_hydrazine^(1/gamma) * ( ((gamma+1)/(gamma-1)) *(1-P_ratio_hydrazine^((1-gamma)/gamma)) )^-.5,P_ratio_hydrazine);

P_ce_ratio_hydrazine = double(P_ce_ratio_hydrazine);

P_ec_ratio_hydrazine = 1 / P_ce_ratio_hydrazine;

CF_hydrazine = (sqrt(   ((2*k_hydrazine^2) / (k_hydrazine-1)) * ( (2/(k_hydrazine+1))^((k_hydrazine+1)/(k_hydrazine-1)) ) * (1 - (P_ec_ratio_hydrazine.^((k_hydrazine-1)/k_hydrazine)))   )) + (P_ec_ratio_hydrazine * Eps);

F_hydrazine = mdot * cstar_hydrazine * CF_hydrazine;

Isp_hydrazine = F_hydrazine / (mdot * g);

% Nitrogen

k_nit = 1.4;
R_nit = 297;
Tc_nit = 300;
cstar_nit = (sqrt(k_nit*R_nit*Tc_nit)) / (k_nit * (sqrt((2/(k_nit+1))^((k_nit+1)/(k_nit-1)))));

syms P_ratio_nit

P_ce_ratio_nit = vpasolve(Eps == ( 2/(k_nit+1))^(1/(k_nit-1)) * P_ratio_nit^(1/k_nit) * ( ((k_nit+1)/(k_nit-1)) *(1-P_ratio_nit^((1-k_nit)/k_nit)) )^-.5,P_ratio_nit);


P_ce_ratio_nit = double(P_ce_ratio_nit);

P_ec_ratio_nit = 1 / P_ce_ratio_nit;

CF_nit = (sqrt(   ((2*k_nit^2) / (k_nit-1)) * ( (2/(k_nit+1))^((k_nit+1)/(k_nit-1)) ) * (1 - (P_ec_ratio_nit.^((k_nit-1)/k_nit)))   )) + (P_ec_ratio_nit * Eps);

F_nit = mdot * cstar_nit * CF_nit;

Isp_nit = F_nit / (mdot * g);

% Hydrogen Peroxide

k_per = 1.25;
R_per = 290;
Tc_per = 1300;

cstar_per = (sqrt(k_per*R_per*Tc_per)) / (k_per * (sqrt((2/(k_per+1))^((k_per+1)/(k_per-1)))));

syms P_ratio_per

P_ce_ratio_per = vpasolve(Eps == ( 2/(k_per+1))^(1/(k_per-1)) * P_ratio_per^(1/k_per) * ( ((k_per+1)/(k_per-1)) *(1-P_ratio_per^((1-k_per)/k_per)) )^-.5,P_ratio_per);


P_ce_ratio_per = double(P_ce_ratio_per);

P_ec_ratio_per = 1 / P_ce_ratio_per;

CF_per = (sqrt(   ((2*k_per^2) / (k_per-1)) * ( (2/(k_per+1))^((k_per+1)/(k_per-1)) ) * (1 - (P_ec_ratio_per.^((k_per-1)/k_per)))   )) + (P_ec_ratio_per * Eps);

F_per = mdot * cstar_per * CF_per;

Isp_per = F_per / (mdot * g);

%%
matlab.codetools.requiredFilesAndProducts("AERO_402_Homework_4_Kulp.m")

