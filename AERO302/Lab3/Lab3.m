%% Lab3
clear all;
close all;
clc

addpath('C:\joshFunctionsMatlab\')
%% Import data
SV1 = importdata("500RPM_6.4deg.mat"); % Scanivalve
SV2 = importdata("500RPM_neg0.8deg.mat");
SV3 = importdata("500RPM_neg5.5deg.mat");
SV4 = importdata("500RPM_neg7.6deg.mat");

DM = [SV1 SV2 SV3 SV4];

clear SV1 SV2 SV3 SV4
rho = 1.204; % kg/m^3

SL1 = importdata("StingMeasuredLoads_20mps_6labdeg.mat");
SL2 = importdata("StingMeasuredLoads_20mps_0labdeg.mat");
SL3 = importdata("StingMeasuredLoads_20mps_neg4labdeg.mat");
SL4 = importdata("StingMeasuredLoads_20mps_neg7labdeg.mat");

SL = [SL1,SL2,SL3,SL4]; % Sting loads each rpm
alpha = [6.4,-.8,-5.5,-7.6]; % raw angle measured
alpha = -alpha; % wing is upside down

LiftMom = importdata('LiftMom.mat'); % Newton.Inches I.G.

c = 4.4482189159;% Newtons/lbs


% Cmat = importdata("FutekCalibrationMatrix-Spring2022.mat"); % Self
Cmat = importdata("X_calibrationFutek_F2022.mat"); % highermath



F = @(dat) dat.F*Cmat - dat.SL.F*Cmat;
Fcal = @(dat) dat.F*Cmat;

% DM force and pressure
for i = 1:length(DM)
    DM(i).SL=SL(i);

    DM(i).P_gauge = mean(DM(i).P); % ?
    DM(i).P = DM(i).P_gauge + DM(i).Amb_Press;
    DM(i).P_tot = DM(i).P(1); % total pressure at test section
    DM(i).P_sta = DM(i).P(1); % static pressure at test section
    DM(i).P_dyn = DM(i).P_tot-DM(i).P_sta; % dynamic pressure
    DM(i).F = DM(i).F(1:3).*c;
    
%     DM(i).res_frc = abs(F(DM(i))); % resultant force
%     DM(i).Fcal = abs(Fcal(DM(i))); % half calibrated?
% 
%     DM(i).F = mean(DM(i).F);
%     DM(i).res_frc = DM(i).res_frc(2:3);
%     DM(i).Fcal = DM(i).Fcal(2:3);

    DM(i).alpha = alpha(i);
    DM(i).LiftMom = LiftMom(i);
end


clear SL SL1 SL2 SL3 SL4 alpha LiftMom

DM = rmfield(DM,'V');
DM = rmfield(DM,"filepath");
DM = rmfield(DM,"t");
DM = rmfield(DM,"Port");

for i = 1:length(DM)
    DM(i).V = sqrt((DM(i).P_dyn*2)/rho); % V at inlet from dynamic pressure
end

for i = 1:length(DM)
    DM(i).SL = rmfield(DM(i).SL,"Amb_Press");
    DM(i).SL = rmfield(DM(i).SL,"V");
    DM(i).SL = rmfield(DM(i).SL,"filepath");
    DM(i).SL = rmfield(DM(i).SL,"Port");
    DM(i).SL = rmfield(DM(i).SL,"Re");
    DM(i).SL = rmfield(DM(i).SL,"t");

    % not sure if P is important
    DM(i).SL = rmfield(DM(i).SL,"P");

    %     DM(i).SL.P_gauge = mean(DM(i).P); % ?
    %     DM(i).SL.P = DM(i).P_gauge + DM(i).Amb_Press;
    %     DM(i).SL.P_tot = DM(i).P(1); % total pressure at test section
    %     DM(i).SL.P_sta = DM(i).P(1); % static pressure at test section
    %     DM(i).SL.P_dyn = DM(i).P_tot-DM(i).P_sta; % dynamic pressure

    DM(i).SL.F = mean(DM(i).SL.F);
    DM(i).SL.F = DM(i).SL.F(1:3).*c;
    p = DM(i).SL.F;
    DM(i).SL.F(1) = p(2);
    DM(i).SL.F(2) = p(1);

    DM(i).res_frc = abs(F(DM(i))); % resultant force
    DM(i).Fcal = abs(Fcal(DM(i))); % half calibrated?

    DM(i).res_frc = DM(i).res_frc(2:3);
    DM(i).Fcal = DM(i).Fcal(2:3);
    DM(i).SL.F = DM(i).SL.F(2:3);
    
end
clear p


c = 0.393701; % cm/in
chord = 4.625; % in
span = 26.5; %in

chord = chord/(c*100); % m
span = span/(c*100); % m
A = chord*span; % m^2

for i = 1:length(DM)
    DM(i).LiftMom = DM(i).LiftMom*c/100; % Newton.meters
end

clear c
AR = span/chord;


%% Lift and Drag resultants

% lift = @(dat) dat.F(2)*Cmat - dat.SL.F(2);
% drag = @(dat) dat.F(1)*Cmat - dat.SL.F(1);



% F = @(dat) Cmat*dat.F' - Cmat*dat.SL.F';
% Fcal = @(dat) Cmat*dat.F';

% for i = 1:length(DM)
%     DM(i).res_frc = abs(F(DM(i))); % resultant force
%     DM(i).Fcal = abs(Fcal(DM(i))); % half calibrated?
% end

clear F
%% Lift and drag coeffs vs angle of attack
Cd = @(dat)dat.res_frc(1)/(.5*rho*A*dat.V^2);
Cl = @(dat)dat.res_frc(2)/(.5*rho*A*dat.V^2);

for i = 1:length(DM)
    DM(i).Cd = Cd(DM(i));
    DM(i).Cl = Cl(DM(i));
    DM(i).LDratio = DM(i).Cl/DM(i).Cd;% lift/drag
end


clear Cd Cl i

%% Coeff Moment
CLM = @(dat) dat.LiftMom/(chord*A*dat.P_dyn);
for i = 1:length(DM)
    DM(i).CLM = CLM(DM(i));
end


%% plots

Fcal = @(dat) dat.F*Cmat;
for i = 1:length(DM)
    Y1(i) = DM(i).Cd;
    Y2(i) = DM(i).Cl;
    Y3(i) = DM(i).LDratio;
    Y8(i) = DM(i).CLM;
%     Y4(i) = DM(i).F(2); % lift frc
%     Y5(i) = DM(i).F(1); % drag frc
% 
%     Y6(i) = DM(i).res_frc(2); % lift frc
%     Y7(i) = DM(i).res_frc(1); % drag frc
%     Y8(i) = DM(i).CLM;
%     Y9(i) = DM(i).Fcal(2);
%     Y10(i) = DM(i).Fcal(1);
    X(i) = DM(i).alpha;
end

%%%%%%%%%%%%% coeffs
figure
plot(X,Y1)
xlabel("$\alpha$ [$^{\circ}$]",'Interpreter','latex')
ylabel("$C_D$",'Interpreter','latex')
% legend("Cd","Cl","CMOM",'Location','best','Interpreter','latex')
title("Coefficient Drag vs $\alpha$",'Interpreter','latex')
grid("on")

figure
plot(X,Y2)
xlabel("$\alpha$ [$^{\circ}$]",'Interpreter','latex')
ylabel("$C_L$",'Interpreter','latex')
% legend("$C_L$",'Location','best','Interpreter','latex')
title("Coefficient Lift vs $\alpha$",'Interpreter','latex')
grid("on")

figure
plot(X,Y8)
xlabel("$\alpha$ [$^{\circ}$]",'Interpreter','latex')
ylabel("$C_M$",'Interpreter','latex')
% legend("$C_M$",'Location','best','Interpreter','latex')
title("Coefficient Moment vs $\alpha$",'Interpreter','latex')
grid("on")

figure
plot(X,Y1,X,Y2,X,Y8)
xlabel("$\alpha$ [$^{\circ}$]",'Interpreter','latex')
ylabel("Coeffs",'Interpreter','latex')
legend("$C_D$","$C_L$","$C_M$",'Location','best','Interpreter','latex')
title("Coefficient Lift, Drag, Moment vs $\alpha$",'Interpreter','latex')
grid("on")



%%%%%%%%%%%%%LIFT/DRAG
figure
plot(X,Y3)
xlabel("$\alpha$ [$^{\circ}$]",'Interpreter','latex')
legend("Lift/Drag ratio",'Location','best','Interpreter','latex')
title("$C_L$/$C_D$ Ratio vs $\alpha$",'Interpreter','latex')
grid("on")


%%%%%%%%%%% forces


Fcal = @(dat) dat.F*Cmat;
for i = 1:length(DM)
    Y1(i) = DM(i).F(2); % lift frc
    Y2(i) = DM(i).F(1); % drag frc
    Y3(i) = DM(i).res_frc(2); % lift frc calibrated-sting
    Y4(i) = DM(i).res_frc(1); % drag frc
    Y5(i) = DM(i).Fcal(2); % calibrated
    Y6(i) = DM(i).Fcal(1);
    X(i) = DM(i).alpha;
end

figure
plot(X,Y1,X,Y2,   X,Y5,X,Y6,    X,Y3,X,Y4)
xlabel("$\alpha$ [$^{\circ}$]",'Interpreter','latex')
ylabel("Force [N]",'Interpreter','latex')
legend("lift force raw","drag force raw","lift force calibrated","drag force calibrated","lift force calibrated minus sting arm","drag force calibrated minus sting arm",'Location','northwest','Interpreter','latex')
title("Lift and Drag vs $\alpha$",'Interpreter','latex')
grid("on")

% figure
% plot(X,Y6,X,Y7)
% xlabel("$\alpha$ [$^{\circ}$]",'Interpreter','latex')
% ylabel("Force [N]",'Interpreter','latex')
% legend("lift force","drag force",'Location','best','Interpreter','latex')
% title("Lift and Drag Raw Values vs $\alpha$",'Interpreter','latex')
% grid("on")









