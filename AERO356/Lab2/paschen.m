clear; close all; clc


%% data read
SSD1 = readtable("AERO 356-06 Arcing Lab Grp 1 Data.xlsx", "NumHeaderLines",2);
P_SSD1 = .5*rmmissing(SSD1.Pressure_Torr_);
V_SSD1 = rmmissing(SSD1.Volts_kV_);

AlD1 = readtable("Arching.xlsx");
P_AlD1 = rmmissing(AlD1.PressureGap_TorrIn_);
V_AlD1 = rmmissing(AlD1.Voltage_V_);

% CONFIGURATION: COPPER ANODE, DISTANCE == 0.5 INCHES
P_CuD1 = .5*[5.9, 3.8, 2.4, 1.2, .28 , .16 ]; % actual values
V_CuD1 = 1e3*[1.000, .700, .700, .600, .621, .784 ]; % voltage in kV when arcing occurred


SSD2 = readtable("Group 4 SS D2.xlsx");
P_SSD2 = SSD2.Pressure_Torr_;
V_SSD2 = SSD2.Voltage_V_;



%% plot setup

labels = {'Stainless Steel, .5" gap', 'Stainless Steel, 1" gap', 'Copper, 1" gap', 'Aluminum, 1" gap'}; %etc

%% plot

hold on
loglog(P_SSD1, V_SSD1*1000, '-')
loglog(P_SSD2, V_SSD2, '-')
loglog(P_CuD1, V_CuD1, '-')
loglog(P_AlD1, V_AlD1, '-')


ax = gca;
ax.XScale = "log";
ax.YScale = "log";

ylim([550 1000])
xlabel("Pressure Gap [Torr*inch]")
ylabel("Breakdown Voltage [V]")
title("Paschen Curves")
legend(labels,"Location","north")