%% 0 
close all;
clear all;
clc

%% Data Processing
addpath('C:\AERO302')
P_R1 = importdata('Group3-200RPM.mat');
P_R2 = importdata('Group3-400RPM.mat');
P_R3 = importdata('Group3-600RPM.mat');
P_R4 = importdata('Group3-800RPM.mat');

vi_R1 = importdata('Group3-200RPM_vi.mat');
vi_R2 = importdata('Group3-400RPM_vi.mat');
vi_R3 = importdata('Group3-600RPM_vi.mat');
vi_R4 = importdata('Group3-800RPM_vi.mat');

vo_R1 = importdata('Group3-200RPM_vo.mat');
vo_R2 = importdata('Group3-400RPM_vo.mat');
vo_R3 = importdata('Group3-600RPM_vo.mat');
vo_R4 = importdata('Group3-800RPM_vo.mat');

Area = importdata('Group3-Area.mat');

T_amb = importdata('Group3-Temp_ambeint.mat');

Power = importdata('Group3-RPM_Freq_Power.mat');

Pres_amb = P_R1.Amb_Press; % Pa

DM = [P_R1,P_R2,P_R3,P_R4];
DM = rmfield(DM,{'filepath','t','Re','Port'});
vi_Mats = [myav(vi_R1),myav(vi_R2),myav(vi_R3),myav(vi_R4)];
vo_Mats = [myav(vo_R1),myav(vo_R2),myav(vo_R3),myav(vo_R4)];

for i = 1:length(DM)
    DM(i).P_gauge = mean(DM(i).P); % ?
    DM(i).P = DM(i).P_gauge + DM(i).Amb_Press;
    DM(i).V = mean(DM(i).V); % ?
%     Data_Mats(i).t = mean(Data_Mats(i).t); % ?

    DM(i).vi = vi_Mats(i); % m/s
    DM(i).vo = vo_Mats(i);

    DM(i).Power = Power(i,3); % Kw
    DM(i).Freq = Power(i,2); % Hz
    DM(i).RPM = Power(i,1); % RPM
    DM(i).T_amb = T_amb; % K
    DM(i).Area = Area; % m^2
    DM(i).P_dyn = DM(i).P(1)-DM(i).P(2);
end


clear Power T_amb vo_Mats vi_Mats Pres_amb vo_R1 vo_R2 vo_R3 vo_R4 vi_R1 vi_R2 vi_R3 vi_R4 P_R1 P_R2 P_R3 P_R4 i

%% Bernoulli Question
% How does your data compare to what Bernoulli would predict if you only had Station 1 to work 
% with? Or station 1 and 2? (really, where does reality start to diverge from Bernoulli in a meaningful 
% way, if it does, and why?).
% bernoulli formula: P/rho + V^2/2 + z*g = const
% neglect z
% rho = P/(R_air*T)
% R_air = R_univesal/(molecularWt)
% R_universal [J/mol/K]
molWt = 28.9647; % g/mol
molWt = molWt/1000; % kg/mol
R_universal = 8.31446261815324; % J/mol/K
R_air = R_universal/molWt; % J/(kg.K)

P = DM(1).Amb_Press;
T= DM(1).T_amb;
% rho = P/(R_air*T); % kg/m^3
rho = 1.204; % kg/m^3

m_dot =@(dat) dat.vi*rho*Area(1); % mass flow rate kg/s
vol_dot =@(dat) dat.Area(1)*dat.vi; % m^3/s
V_fromA = @(dat) (dat.Area.^-1).*dat.vol_dot; % calculates Velocity (m/s) from an Areas vector and vol_dot

for i = 1:length(DM)
    DM(i).rho = mean(rho); % kg / m^3
    DM(i).m_dot = m_dot(DM(i)); % kg / s
    DM(i).vol_dot = vol_dot(DM(i)); % m^3 / s
    DM(i).V_fromA = V_fromA(DM(i));

end

% f = figure;
% for i = 1:length(DM)
%     figure
%     nexttile
%     X=categorical({'1','2','3','4','5'});
%     Y = DM(i).V_fromA';
%     bar(X,Y)
%     title(string(DM(i).RPM)+" RPM")
% %     set(gca, 'Ydir','reverse')
%     ylim([0 40])
%     ylabel("Velocity [m/s]")
%     xlabel("Static Pressure Station")
% end
% set(f,'NumberTitle','off','Name',"Velocity as Calculated by Volume Flow Rate and Area")

clear i P R_air R_universal molWt m_dot f


cfun = @(dat) dat.P(2)/dat.rho + (dat.vi)^2/2; % this is the constant value that bernoulli adds to: P/rho + V^2/2 = const
P_perdict = @(dat) ((-(dat.V_fromA).^2)/2+cfun(dat))*dat.rho; % bernoulli

for i = 1:length(DM)
    DM(i).P_perdict = P_perdict(DM(i));
    DM(i).P_descrepency = DM(i).P((2:6))-DM(i).P_perdict;
    DM(i).P_perdict_gauge = DM(i).P_perdict-DM(i).Amb_Press;
end

% f = figure;

% for i = 1:length(DM)
%     figure
%     t = tiledlayout(2,2);
%     hold on
%     nexttile
%     X=categorical({'1','2','3','4','5'});
%     Y = [(DM(i).P_gauge((2:6))'),(DM(i).P_perdict_gauge'),(DM(i).P_descrepency')];
%     bar(X,Y)
%     title(string(DM(i).RPM)+" RPM")
%     set(gca, 'Ydir','reverse')
%     ylim([-1000 250])
%     ylabel("Pressure [Pa]")
%     xlabel("Static Pressure Station")
%     legend("Measured P" ,"Bernoulli P","Difference",'Location','bestoutside');
% end



% legend("Measured P" ,"Bernoulli P","Difference",'Location','bestoutside');
% set(f,'NumberTitle','off','Name',"Comparison of Measured and Predicted Gauge Pressure Drop [Pa]")


%%  Coeff Pressure
clear cfun f t i l P_perdict V_fromA vol_dot X Y T rho

rho = DM(1).rho(1);
Pamb = DM(1).Amb_Press(1);
CP=@(dat) 2*(dat.P(2:6)-Pamb)./(rho*dat.V_fromA.^2); % formula for coeff of pressure
for i = 1:length(DM)
    DM(i).CP = CP(DM(i));
end

% % f = figure;
% for i = 1:length(DM)
%     figure
%     nexttile
%     X=categorical({'1','2','3','4','5'});
%     Y = DM(i).CP';
%     bar(X,Y)
%     title(string(DM(i).RPM)+" RPM")
%     set(gca, 'Ydir','reverse')
%     ylim([-20 0])
%     ylabel("Coefficient of Pressure")
%     xlabel("Static Pressure Station")
% end
% set(f,'NumberTitle','off','Name',"Coefficient or Pressure at each station")

%% exhaust speed


%% Functions
function [av] = myav(vec)
    vec = sort(vec);
    vec(end) = [];
    vec(1) = [];
    av = mean(vec);
end

%%
