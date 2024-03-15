%% 
clear all
close all
clc

%% morhs circle
sigy = 10;
sigx = 0;
c = (sigy+sigx)/2;
tmax = (sigy-sigx)/2;

figure
x = [sigx, sigy, c]
y = [0,0,tmax]

%% data
close all

D = readmatrix('right/Test1/Test1.steps.trends.csv');
n = 2936
load = mean(D(1:n,6)-D(1:n,7))
K = 2.45
A = (1)*(1/8) * 0.00064516 % m^2
sigAve = load/A;
sigMax = sigAve*K
sigMax = sigMax/1e3 % Mpa

a = 7.10 
w = 25.4
Af = A*(a/w)
Scrit = load/A;
Scrit = Scrit/1e3

Kic = 29e6
acrit = ((Kic/(sigAve))^2 )/pi;
acrit = acrit * 1e3
scatter(n,sigMax)
ylim([0,500])
xscale("log")

% D = D(1000:2400,:);
% 
% E = 69e9; %pa = 69Gpa
% % sig = E * epsilon
% len = 8; % cm
% epsilon = D(:,4) / len*10;
% sig = E * epsilon;
% 
% figure
% hold on
% plot(D(:,1),sig)
% xlabel("cycles")
% ylabel("sigma")
% xscale("log")
% 
% 
% figure
% hold on
% plot(D(:,1),D(:,4))
% xlabel("cycles")
% ylabel("Position:Maximum (mm)")
% xscale("log")
% 
% header = ["Total Cycles","	Elapsed Cycles","Step","Position:Maximum (mm)","Position:Minimum (mm)","Load:Maximum (kN)","Load:Minimum (kN)","Load:Amplitude (kN)","Load:Mean of Peaks (kN)"];
% for i = 4:9
%     figure
%     hold on
%     plot(D(:,1),D(:,i))
%     xlabel("cycles")
%     ylabel(header(i))
%     xscale("log")
% end
