%% Joshua Oates - Estes Rocket Lab
clear all
close all
clc


%% import
D12weight = importdata("D12weight.mat");
C6weight = importdata("C6weight.mat");

D12fire = importdata("D12fire.mat");
C6fire = importdata("C6fire.mat");

%% calculate and remove the unloaded average
dt = .2; %s, time step of measurements
avet = 5; %s, time that we average unloaded data over
avepts = avet/dt;
D12unload = mean(D12fire(1:avepts,2));
C6unload = mean(C6fire(1:avepts,2));

% remove unloaded average from data, remove unused vars
D12fire(:,2) = D12fire(:,2)-D12unload;
C6fire(:,2) = C6fire(:,2)-C6unload;
clear D12unload C6unload avepts avet

% remove data surrounding firing
C6tspan = [6.6,9.8]; % seconds
ptspan = C6tspan./dt;
ptspan = round(ptspan);
C6fire = C6fire(ptspan(1):ptspan(2),:);

D12tspan = [7.2,9.8]; % seconds
ptspan = D12tspan./dt;
ptspan = round(ptspan);
D12fire = D12fire(ptspan(1):ptspan(2),:);
clear ptspan

% renumber the seconds colunm
lin = 0:dt:length(C6fire)*dt-dt;
C6fire(:,1) = lin';

lin = 0:dt:length(D12fire)*dt-dt;
D12fire(:,1) = lin';
clear lin

% lbs to N
N_lbs = 4.4482189159; % = 1lbs
D12fire(:,2) = D12fire(:,2).*N_lbs;
C6fire(:,2) = C6fire(:,2).*N_lbs;
clear N_lbs

%% Ploting

figure
plot(C6fire(:,1),C6fire(:,2))
ylim([0,15])
xlim([0,3])
title("C6")
xlabel("time [s]")
ylabel("thrust [N]")


figure
plot(D12fire(:,1),D12fire(:,2))
ylim([0,30])
xlim([0,3])
title("D12")
xlabel("time [s]")
ylabel("thrust [N]")

%% calcs avg thrust total impulse

C6ave = mean(C6fire(:,2))
D12ave = mean(D12fire(:,2))

C6I = sum(C6fire(:,2)*dt)
D12I= sum(D12fire(:,2)*dt)

C6dm = C6weight(1)-C6weight(2)
D12dm = D12weight(1)-D12weight(2)

C6mdot = C6dm/(C6tspan(2)-C6tspan(1))
D12mdot = D12dm/(D12tspan(2)-D12tspan(1))

