%% AERO 302 - Lab 2 scratch

%% 0 
clear all; close all; clc

%% Imports
SV1 = importdata("cylindertest1.mat"); % Scanivalve
SV2 = importdata("cylindertest2.mat"); 
SV3 = importdata("cylindertest3.mat");
SV4 = importdata("cylindertest4.mat");
SV5 = importdata("cylindertest5.mat");
SV6 = importdata("cylindertest6.mat");
V_sc = importdata("V_sc.mat");
DM = [SV1 SV2 SV3 SV4 SV5 SV6];

clear SV1 SV2 SV3 SV4 SV5 SV6
rho = 1.204; % kg/m^3
% from lab view the first port is tot second is stat, 3-26 is cylinder
rho = 1.204; % kg/m^3
for i = 1:length(DM)

    DM(i).std = std(DM(i).P); % std dev
        DM(i).covar = DM(i).std./mean(abs(DM(i).P(:)));
    DM(i).V_sc = V_sc(i); % speed controller
    DM(i).V_sv = DM(i).V; % scanivalve
    DM(i).V_sv = mean(DM(i).V_sv);
    DM(i).P_gauge = mean(DM(i).P); % ?
    DM(i).P = DM(i).P_gauge + DM(i).Amb_Press;
    DM(i).P_tot = DM(i).P(1); % total pressure at test section
    DM(i).P_sta = DM(i).P(2); % static pressure at test section
    DM(i).P_dyn = DM(i).P_tot-DM(i).P_sta; % dynamic pressure
    DM(i).P_gauge = DM(i).P_gauge(3:end);
    DM(i).P = DM(i).P(3:end);

end
RPM = ['500 RPM sym','600 RPM sym','300 RPM sym','500 RPM trip','600 RPM trip','300 RPM trip'];
%RPM = string(RPM);

DM(1).RPM = '500 sym';
DM(2).RPM = '600 sym';
DM(3).RPM = '300 sym';
DM(4).RPM = '500 trip';
DM(5).RPM = '600 trip';
DM(6).RPM = '300 trip';


clear RPM
DM = rmfield(DM,'V');
for i = 1:length(DM)
    DM(i).V = sqrt((DM(i).P_dyn*2)/rho); % V at inlet from dynamic pressure
end
DM = rmfield(DM,"filepath");
DM = rmfield(DM,"t");
clear V_sc
% 
% figure
% Y = DM(5).P_gauge;
% bar(Y)
% figure
% Y = DM(2).P_gauge;
% bar(Y)
%% setting up port angles
n = 24;
ang_inc = (2*pi)/24;
prt_ang=zeros(1,n); % port angle
prt_unt=zeros(2,n); % port unit vector
for i = 1:n
    prt_ang(i) = ang_inc*(i-1); % here index 1 goes with port 24 since port 24 is nominally on the x axis
end 
prt_ang = prt_ang+ang_inc; % now the angles can map directly to the 
prt_ang(n) = 0;
for i = 1:n
    prt_unt(:,i) = -[cos(prt_ang(i));sin(prt_ang(i))]; % inward facing vectors that will be in the direction of force
end
dia = 155; % diameter in mm
dia = dia/1000;
cir = pi*dia; % circumfrenc in m
len = cir/24; % length of each section
A = 1*len; % area of each panel m^2, asumme cylinder is 1m long
clear dia cir len
for i = 1:length(DM)
    DM(i).prt_frc=zeros(2,n); % force vector
    P1 = DM(i).P;
    P2 = P1;
    P2(23) = (P2(22)+P2(24))/2; % 23 is an outlyer
    DM(i).P2 = P2;
    
    for j =1:n
%         DM(i).prt_frc1(:,j) = prt_unt(:,j).*P1(j)*A; % force vector acting per unit area on each port
        DM(i).prt_frc2(:,j) = prt_unt(:,j).*P2(j)*A;
    end
%     DM(i).res_frc1 = [sum(DM(i).prt_frc1(1,:));sum(DM(i).prt_frc1(2,:))];
    DM(i).res_frc2 = [sum(DM(i).prt_frc2(1,:));sum(DM(i).prt_frc2(2,:))];
end
clear A P1 P2 i j ang_inc %prt_unt prt_ang
%% Cp vs theta
for i = 1:length(DM)
    P_inf = DM(i).P_sta; % freestream absolute static pressure, an averaged value from station 2 
    V = DM(i).V; % freestream Velocity
    CP = @(P) (P-P_inf)/(.5*rho*V^2); % definition of CP(P) where P is absolute pressure measured at a given port
    for j = 1:n
        DM(i).CP(j) = CP(DM(i).P(j)); % calculate CP for all ports
    end        
end
prt_ang(n) = 2*pi;

%%%% reorganize DM
sym300 = DM(3);
sym500 = DM(1);
sym600 = DM(2);
trip300 = DM(6);
trip500 = DM(4);
trip600 = DM(5);

DM(1) = sym300;
DM(2) = sym500;
DM(3) = sym600;
DM(4) = trip300;
DM(5) = trip500;
DM(6) = trip600;

%%% Cp theoretical

CpT = @(prt_ang) 1-4*(sin(prt_ang)).^2;

%%%%%% figure
figure
hold on
grid on
    % sym
    plot(prt_ang*(180/pi),DM(1).CP,'o-')
    plot(prt_ang*(180/pi),DM(2).CP,'o-')
    plot(prt_ang*(180/pi),DM(3).CP,'o-')
    % trip
    plot(prt_ang*(180/pi),DM(4).CP,'*-')
    plot(prt_ang*(180/pi),DM(5).CP,'*-')
    plot(prt_ang*(180/pi),DM(6).CP,'*-')

    plot(prt_ang*(180/pi), CpT(prt_ang), 's-')


title('Measured and Theoretical Coefficient of Pressure versus Port Angle','Interpreter','latex')
xlabel('$\theta$ (deg)','Interpreter','latex')
ylabel('$C_{p}$','Interpreter','latex')
legend(DM.RPM,'$C_{P}$ Theoretical','Interpreter','latex','Location','Northwest')
xlim([0 360])


%%%%%%%%%%%redo with P2
for i = 1:length(DM)
    P_inf = DM(i).P_sta;
    V = DM(i).V;
    CP = @(P) (P-P_inf)/(.5*rho*V^2);
    for j = 1:n
        DM(i).CP(j) = CP(DM(i).P2(j));
    end        
end
prt_ang(n) = 2*pi;

% figure
% hold on
% grid on
% for i = 1:length(DM)
%     plot(prt_ang,DM(i).CP)
% end
% title('Measured Coefficient of Pressure versus Port Angle','Interpreter','latex')
% xlabel('$\theta$ (rad)','Interpreter','latex')
% ylabel('$C_{p}$','Interpreter','latex')
% legend(DM.RPM,'Location','Northwest')


%%%%%% STD DEV
% plot coefficients of variation
figure()

grid on
hold on
for i = 1:length(DM)
    plot(DM(i).covar, 'o', 'LineWidth',2)
end
title('Coefficient of Variation at Each Pressure Port','Interpreter','latex')
xlabel('Pressure Port Number','Interpreter','latex')
ylabel('Coefficient of Variation','Interpreter','latex')
legend(DM.RPM,'Interpreter','latex','Location','northwest')


%%%%%%% RE num
figure()
grid on
hold on

barx = categorical({'∆Re 300 RPM','∆Re 500 RPM','∆Re 600 RPM'});
bary = [DM(4).Re-DM(1).Re, DM(5).Re-DM(2).Re, DM(6).Re-DM(3).Re];
bar(barx, bary)

title('RPM vs. Reynolds Number with and without Trip Strip','Interpreter','latex')
xlabel('RPM','Interpreter','latex')
ylabel('Re','Interpreter','latex')



%% Functions
function [av] = myav(vec)
    vec = sort(vec);
    vec(end) = [];
    vec(1) = [];
    av = mean(vec);
end


