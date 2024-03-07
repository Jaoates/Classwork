clear all;
close all;
clc
addpath('C:\joshFunctionsMatlab\')

h = 10*1000;
[T,P,rho,mu] = joshStdAtm(h)

gamma = 1.4;
Rgas = 287.85; % J kg K
a = sqrt(gamma*Rgas*T)

M = .5 % = U/a
U = M*a

[TRatio,PRatio,rhoRatio]=joshIsentropicToolbox(M)

T1 = T/TRatio
P1 = P/PRatio
rho1 = rho/rhoRatio

muRef = 1.716e-5;
TRef = 273.15;
S = 110.4;


mufun = @(T) muRef*(((T+S).^-1).*(TRef+S)).*(T/TRef).^1.5;

mu = mufun(T1)
% mu1 = mufun(T1)

Re = rho*U*2/mu

cp = (2.6431e4-2.2282e4)/(.5*rho*U^2)

%% 2
clear all
[T,P,rho,mu] = joshStdAtm(.1)

C = sqrt(P/rho)
Rgas = 287
T = 300
gamma = 1.4
a = sqrt(gamma*Rgas*T)
