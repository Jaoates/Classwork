clear all
close all
clc

addpath('C:\AERO303\compressibleFlowRelations\')
%% 
gam = 1.4;
R = .287;%kJ/(kg*K)
R = R*1000;%J/(kg*K)

cp = 1.005;%kJ/(kg*K)
cp = cp*1000;%J/(kg*K)

T1 = 269;
P1 = 70000;
M1 = 1.8;
n1 = joshComp(M1,2,'M')

P2 = n1.P2_P1*P1
T2 = n1.T2_T1*T1

ds = cp*log(n1.T2_T1)-R*log(n1.P2_P1)

%%
gam = 1.4;
R = .287;%kJ/(kg*K)
R = R*1000;%J/(kg*K)

cp = 1.005;%kJ/(kg*K)
cp = cp*1000;%J/(kg*K)

M1=.9
P1 = 1;%atm
T1 = 600;%K
q = -578078;%J/kg

isen1 = joshComp(M1,1,'M')
P01 = isen1.P0_P*P1;
T01 = isen1.T0_T*T1;

T02 = (q+cp*T01)/cp

heat1 = joshComp(M1,1,'M');




