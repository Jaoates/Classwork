%% Props hw5 - Joshua Oates
clear all
close all
clc

%% Problem 1

disp("Problem 1")
epsilon = 60 
Mpay = 20000 
dv = 500 
t = 10*60 
MM = [1.00784,15.999]*2 
CoeffR = [2,1] 
MMR = MM.*CoeffR 
r = [MMR(2)/MMR(1),3.8,6] 
MMC = [16,9.5,13.2] 
R = 8.31432 *1000 
R = R./MMC 
% gam = Cp/Cv
gam = [1.2,1.225,1.2] 
Tc = [3500,2800,3400] 
% ISPvac = [386,435,414] 
syms Pce
eqn1 = epsilon == (2./(gam+1)).^(1./(gam-1)).*Pce.^(1./gam).*(((gam+1)./(gam-1)).*(1-Pce.^((1-gam)./gam))).^-.5 % definition expansion ratio lecture3
Pce = [vpasolve(eqn1(1),Pce),vpasolve(eqn1(2),Pce),vpasolve(eqn1(3),Pce)] 
Pce = double(Pce) 
Pec = 1./Pce 
Cf = sqrt((2*gam.^2./(gam-1)) .* ((2./(gam+1)).^ ((gam+1)./(gam-1)) .* (1-Pec.^((gam-1)./gam)))) + Pec.*epsilon 
g = 9.81 
Cs = sqrt(gam.*R.*Tc)./(gam.*sqrt((2./(gam+1)).^((gam+1)./(gam-1)))) % C star from lecture7
Isp = Cs.*Cf/g 
syms M0
eqn2 = dv == -Isp.*g.*log(Mpay/M0) 
M0 = [vpasolve(eqn2(1),M0),vpasolve(eqn2(2),M0),vpasolve(eqn2(3),M0)] 
M0 = double(M0) 
Mprop = M0-Mpay 
syms MO2 MH2
eqn3 = [Mprop == MO2+MH2;r == MO2/MH2] 
sol = [solve(eqn3(:,1)),solve(eqn3(:,2)),solve(eqn3(:,3))] 
[MO2(1),MO2(2),MO2(3)] = sol.MO2 
[MH2(1),MH2(2),MH2(3)] = sol.MH2 
MO2 = double(MO2) 
MH2 = double(MH2) 
rho = [70.85,1.141*1000] 
V = [MH2;MO2]./rho' 
VH2 = V(1,:) 
VO2 = V(2,:) 
Cfp =sqrt(((2*gam.^2)./(gam-1)).*( (2./(gam+1)).^((gam+1)./(gam-1)) )) 
IspIdeal = Cs.*Cfp./g 
It = Isp.*Mprop.*g 
mdot = Mprop./t 
F = Isp.*mdot.*g 

%% Problem 2

disp("Problem 2")
clear 
g = 32.1741 
Pe = 14.6959 
Pc = 1000 
gam = 1.26 
Tc = 2700 
br = .1 
Cs = 4000 
rho = .056 
MM = 22 
F = 2000 
t = 10 
Tamb = 70 
IWr = 143 
eta = .98 

Pce = Pc/Pe 
Pec = 1/Pce 
Cf = sqrt((2*gam.^2./(gam-1)) .* ((2./(gam+1)).^ ((gam+1)./(gam-1)) .* (1-Pec.^((gam-1)./gam))))  
Cf = Cf*eta 
At = F/(Cf*Pc) 
epsilon = (2./(gam+1)).^(1./(gam-1)).*Pce.^(1./gam).*(((gam+1)./(gam-1)).*(1-Pce.^((1-gam)./gam))).^-.5 % definition expansion ratio lecture3
Ae = At*epsilon 
Isp = Cs*Cf/g 
mdot = F/Isp 
Ab = Pc*At/((Cs*12)*br*rho) 
Mprop = mdot*t*1.04 
It = Isp*Mprop 
M0 = It/IWr 
Ms = M0-Mprop 

%% Problem 3
disp("Problem 3")
clear
g = 9.81 
F = 7500 
t = 100 
z = 10 
P0 = 26.5e3 
Pc = 10e6 
% item d is clean 
a = .41 
n = .39 
Cs = 1511 
rho = 1680 
Tc = 3288 
gam = 1.2 
Pe = P0 
Pec = Pe/Pc 
Pce = 1./Pec 
Cf = sqrt((2*gam.^2./(gam-1)) .* ((2./(gam+1)).^ ((gam+1)./(gam-1)) .* (1-Pec.^((gam-1)./gam)))) 
At = F/(Cf*Pc) 
Dt = 2*sqrt(At/pi) 
epsilon = (2./(gam+1)).^(1./(gam-1)).*Pce.^(1./gam).*(((gam+1)./(gam-1)).*(1-Pce.^((1-gam)./gam))).^-.5 
Ae = At*epsilon 
De = 2*sqrt(Ae/pi) 
rdot = a*(Pc/1e6)^n 
rdot = rdot/100 
Ab=Pc*At/(Cs*rdot*rho) 
mdot = Ab*rho*rdot 
thickness = rdot*t 
Mprop= mdot*t 
Isp = F/(mdot*g)



%% Problem 4
disp("Problem 4")
close all
d = 20 
d = d/100 
l = 50 
l = l/100 
dp = 4 
dp = dp/100 
rp = dp/2 
rho = 920 
a = .0018 
n = 0.3 
t = linspace(0,100,1000) ;
mdotox = 4 

rp = (a*(2*n+1) * t*(mdotox/pi).^n + rp^(2*n+1)).^(1/(2*n+1)); 
rp = rp(rp<=d/2) ;
t = t(1:length(rp)); 

Ap = pi.*rp.^2 ;
rdot = a.*(mdotox./Ap).^n ;
mdotfuel = rho.*pi.*rp.*rdot; 
figure
plot(t,rdot)
title("Regression rate over time")
xlabel("time [s]")
ylabel("Regression rate [m/s]")
figure
hold on
plot(t,mdotfuel)
plot([t(1),t(end)],[1,1]*mdotox)
plot(t,mdotfuel+mdotox)
legend(["fuel","ox","total"])
ylabel("mass flow rate [kg/s]")
xlabel("time [s]")
title("Mass Flow Rate vs Time")
