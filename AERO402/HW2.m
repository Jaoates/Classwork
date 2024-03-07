%% HW2 - Prop
clear all
close all
clc

%% Problem 2
disp("Problem 2")
g = 9.81;
Arat = [1,1,2,5,10,100];
Cf = [.6636,.6636,1.2157,1.4872,1.6253,1.8976];
Isp = [1562.2,1562.2,2861.9,3501,3826,4467.1];
Isp = Isp./g;

figure
plot(Arat,Cf)
title("Cf and Area Ratio")
xlabel("Epsilon")
ylabel("Cf")

figure
plot(Arat,Isp)
title("Isp and Area Ratio")
xlabel("Epsilon")
ylabel("Isp [s]")

%% Problem 3
disp("Problem 3")
gam = 1.3;
Pc = 20; 
Arat = 6;

Pe = 20/50.5;

Psea = 1; 
P25 = 0.025;  % from std atm model

Pec = Pe/Pc;
Po = P25;
myVar = ( Pe-Po )/Pc *Arat;
Cf(1) = sqrt((2*gam^2/(gam-1)) * ((2/(gam+1))^ ((gam+1)/(gam-1)) * (1-Pec^((gam-1)/gam)))) + myVar;
 
Po = Psea;
myVar = ( Pe-Po )/Pc *Arat;
Cf(2) = sqrt((2*gam^2/(gam-1)) * ((2/(gam+1))^ ((gam+1)/(gam-1)) * (1-Pec^((gam-1)/gam)))) + myVar;

variation = ((Cf(1)-Cf(2))/Cf(2))


%% Problem 4
disp("Problem 4")
% part a
clear
load("IntegralOfCp.mat")
labels = ["CO","CO2","H2","H2O","N2"];

myIntCO = @(Tc) -(ppval(P{1},298)-ppval(P{1},Tc));
myIntCO2 = @(Tc) -(ppval(P{2},298)-ppval(P{2},Tc));
myIntH2 = @(Tc) -(ppval(P{3},298)-ppval(P{3},Tc));
myIntH2O = @(Tc) -(ppval(P{4},298)-ppval(P{4},Tc));
myIntN2 = @(Tc) -(ppval(P{5},298)-ppval(P{5},Tc));

myVecFun = {myIntCO,myIntCO2,myIntH2,myIntH2O,myIntN2};

np = [.836 .164 .664 .836 .5]; % mol
dHp = [-110.53 -393.522 0 -241.826 0]; %kJ/mol
dHr = -113.1;


dHrxn = sum(np.*dHp) - dHr

Tc = 200:1:6000;
myDiff =[];
for i = 1:length(Tc)
    val(i) = 0;
    for j = 1:5
        val(i) = val(i) + np(j)*myVecFun{j}(Tc(i));
    end
end

myDiff = dHrxn+val;
[m,i] = min(abs(myDiff));
Tc = Tc(i);


Tc = Tc-1:.001:Tc+1;
myDiff =[];
for i = 1:length(Tc)
    val(i) = 0;
    for j = 1:5
        val(i) = val(i) + np(j)*myVecFun{j}(Tc(i));
    end
end

myDiff = dHrxn+val;
[m,i] = min(abs(myDiff));


Tc = Tc(i)

% part b

MM = [28.01,44.01,2.016,18.01528,28.02]; % g/mol
MM = MM./1000; %kg/mol

massMixture = sum(MM.*np)/sum(np) % kg/mol

% part c
val = 0;
for j = 1:5
    Cp(j) = (-myVecFun{j}(Tc))/(298-Tc); % kJ/(mol K)
end

CpMixture = sum(np.*Cp)/sum(np); % kJ/(mol K)
CpMixture = CpMixture/massMixture % kJ/(kg K)


% part d

Ru = 8.314; % J/mol-K
Rgas = Ru/massMixture; % (J/mol-K)/(kg/mol) = J/kg-K
Rgas = Rgas/1000 % kJ/kg-K
% part e

CvMixture = CpMixture - Rgas

gamMixture = CpMixture/CvMixture

% extra part
Pce = 69;
Pec = 1/Pce;
gam = gamMixture;

g = 9.81;

Cf = sqrt((2*gam^2/(gam-1)) * ((2/(gam+1))^ ((gam+1)/(gam-1)) * (1-Pec^((gam-1)/gam))));
Cs = (sqrt(gam*Rgas*1000*Tc)) / (gam*sqrt((2/(gam+1))^((gam+1)/(gam-1))));
Isp=Cs*Cf/g

%% Problem 5
disp("Problem 5")
clear all
syms  Prat 

Pco = [2,10,100,1000,inf];
Arat = logspace(0,2);
gam = 1.2;

figure
hold on 

for j = 1:length(Pco)
    for i = 1:length(Arat)
        Pce=vpasolve(Arat(i) == ( 2/(gam+1))^(1/(gam-1)) * Prat^(1/gam) * ( ((gam+1)/(gam-1)) *(1-Prat^((1-gam)/gam)) )^-.5,Prat);
        Poc = 1/Pco(j);
        Pec = 1/Pce;
        myVar = ( Pec - Poc ) *Arat(i);
        Cf(i) = sqrt((2*gam^2/(gam-1)) * ((2/(gam+1))^ ((gam+1)/(gam-1)) * (1-Pec^((gam-1)/gam)))) + myVar;
    end
    [M,I]=max(Cf);
    plot(Arat,Cf)
    scat(:,j) = [Arat(I),Cf(I)];
end


scatter(scat(1,:),scat(2,:))
set(gca,"XScale","log")
legend([("P1/P3 = "+["2","10","100","1000","inf"]) , ["Best Cf"] ],"Location","best")
ylabel("Cf")
xlabel("Epsilon")
ylim([0,2])
