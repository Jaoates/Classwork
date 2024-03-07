%% housekeeping
clear all
close all
clc

%% Load data
d = readmatrix("StrainData2024_02_13.TXT")
d = d(:,end-2:end-1)
d = rmmissing(d)
d = d(4:end,:)

closedtill = 12
openedat = 14

epsilon(1,:) = mean(d(1:closedtill,:))
epsilon(2,:) = mean(d(openedat:end,:))

epsilon = epsilon(2,:)-epsilon(1,:)
strainRatio = epsilon(1)/epsilon(2)

epsilon = epsilon*1e-6;
E = 69e9; %Gpa

sigma = epsilon*E*1e-6
stressRatio = sigma(1)/sigma(2)

vonMises = sqrt((sigma(1)-sigma(2))^2 + (sigma(2)-0)^2 + (0-sigma(1))^2)
FS = 170/vonMises

t = .109e3
D = 65.64e3
nu = .325
P = (4*t*E/D) * [epsilon(1)/(2-nu),epsilon(2)/(1-2*nu)];
P = P*1e-3
P = mean(P)

