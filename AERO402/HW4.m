%% HW 4 - Props
clear all
close all
clc

%% Problem 1
disp("Problem 1")
disp("--single stage--")
Ve = 3048;
Mpay = 1000;
M0 = 15000;
Ms = 2000;
Mprop = M0 - Ms - Mpay;
disp("M0: "+string(M0)+"kg")
disp("Mpay: "+string(Mpay)+"kg")
disp("Ms: "+string(Ms)+"kg")

epsilon = Ms/(Mprop+Ms);
disp("Structural Coeff (epsilon): "+string(epsilon))
lamda = Mpay/(Mprop+Ms);
disp("Payload Ratio (lamda): "+string(lamda))
R = M0/(Ms+Mpay);
disp("Mass Ratio (R): "+string(R))
c = Ve;
Dv = -c*log((Mpay+Ms)/M0);
disp("Delta V: "+string(Dv)+"m/s")
finalV = Dv;
disp("Final Velocity: "+string(finalV)+"m/s")

%Staged part
syms M01 M02 Ms1 Ms2 Mpay1 Mpay2 Mprop1 Mprop2
M01 = M0;
Mpay2 = Mpay;
M02 = Mpay2 + Mprop2 + Ms2;
Mpay1 = M02;
eqn1 = Ms1/(Mpay1+Ms1) == Ms2/(Mpay2+Ms2);
eqn2 = Mpay1/(Mprop1+Ms1) == Mpay2/(Mprop2+Ms2);
eqn3 = M01/(Ms1+Mpay1) == M02/(Ms2+Mpay2);
eqn4 = Ms1+Ms2 == Ms;
eqns = [eqn1,eqn2,eqn3,eqn4];

assume(Ms2>0)
sol = solve(eqns);
% Ms1 = sol.Ms1;
% Ms2 = sol.Ms2;
% Mprop1 = sol.Mprop1;
% Mprop2 = sol.Mprop2;

Dv2 = -c*log((Mpay2+Ms2)/M02);
Dv1 = -c*log((Mpay1+Ms1)/M01);
things = [Dv2,Dv1,M01,M02,Ms1,Ms2,Mpay1,Mpay2,Mprop1,Mprop2];
things = subs(things,Mprop1,sol.Mprop1);
things = subs(things,Mprop2,sol.Mprop2);
things = subs(things,Ms1,sol.Ms1);
things = subs(things,Ms2,sol.Ms2);
Dv1 = double(things(1));
Dv2 = double(things(2));
M01= double(things(3));
M02= double(things(4));
Ms1= double(things(5));
Ms2= double(things(6));
Mpay1= double(things(7));
Mpay2= double(things(8));
Mprop1= double(things(9));
Mprop2= double(things(10));

%%%%
disp("--first stage--")
disp("M0: "+string(M01)+"kg")
disp("Mpay: "+string(Mpay1)+"kg")
disp("Ms: "+string(Ms1)+"kg")

epsilon = Ms1/(Mprop1+Ms1);
disp("Structural Coeff (epsilon): "+string(epsilon))
lamda = Mpay1/(Mprop1+Ms1);
disp("Payload Ratio (lamda): "+string(lamda))
R = M01/(Ms1+Mpay1);
disp("Mass Ratio (R): "+string(R))
things = -c*log((Mpay1+Ms1)/M01);
disp("Delta V: "+string(Dv1)+"m/s")
%%%%
disp("--second stage--")
disp("M0: "+string(M02)+"kg")
disp("Mpay: "+string(Mpay2)+"kg")
disp("Ms: "+string(Ms2)+"kg")

epsilon = Ms2/(Mprop2+Ms2);
disp("Structural Coeff (epsilon): "+string(epsilon))
lamda = Mpay2/(Mprop2+Ms2);
disp("Payload Ratio (lamda): "+string(lamda))
R = M02/(Ms2+Mpay2);
disp("Mass Ratio (R): "+string(R))
things = -c*log((Mpay2+Ms2)/M02);
disp("Delta V: "+string(Dv2)+"m/s")

disp("--staged--")
disp("Final Velocity: "+string(Dv1+Dv2)+"m/s")

%% Problem 2
clear
% close all
disp("Problem 2")
F = 40;
It = 40000;
gam = [1.4,1.2];
Tc = [298,1650];
R = [297,260];
labels=["N2","MMH"];
g = 9.8;
Cs = sqrt(gam.*R.*Tc)./(gam.*sqrt((2./(gam+1)).^((gam+1)./(gam-1))));% C star from lecture7

% ideal Cf vs Pec
Pce = 1:3000;
Pec = 1./Pce;%Pe/Pc
Cf =@(gam) sqrt(... % Coeffecient of thrust from lecture7 for Ideal nozzle
    ((2*gam.^2)./(gam-1)).*...
    ( (2./(gam+1)).^((gam+1)./(gam-1)) ).*...  
    (1-Pec.^((gam-1)./gam)) );

Cf = [Cf(gam(1))',Cf(gam(2))'];
figure
plot(Pce,Cf(:,1),Pce,Cf(:,2))
legend(labels)
title("Cf vs Pc/Pe - Ideal")
xlabel("Pc/Pe")
ylabel("Cf")

% non ideal Cf vs epsilon
epsilon = (2./(gam+1)).^(1./(gam-1)).*Pce'.^(1./gam).*(((gam+1)./(gam-1)).*(1-Pce'.^((1-gam)./gam))).^-.5;
mdot = F./(Cs.*Cf);
Isp = F./(mdot.*g);
figure
plot(epsilon(:,1),Isp(:,1),epsilon(:,2),Isp(:,2))
legend(labels)
title("Isp vs Epsilon")
xlabel("Epsilon")
ylabel("Isp [s]")

%set epsilon to the one found in class
epsilonReal = 25;
[m,i] = min(abs(epsilon-epsilonReal));
Isp = [Isp(i(1),1),Isp(i(2),2)];
epsilon = epsilonReal;
Cf = Isp.*g./Cs;
mdot = F./(Cs.*Cf);
Pc = 0:1e5:2e7;
At = mdot.*Cs./(Pc');
Ae = epsilon.*At;
% At = pi*(Dt/2)^2;
Dt = sqrt(At./pi).*2;
De = sqrt(Ae./pi).*2;
figure
plot(Pc,Dt(:,1),Pc,De(:,1),Pc,Dt(:,2),Pc,De(:,2))
legend("N2 De","N2 Dt","MMH De","MMH Dt")
xlabel("Pc [pa]")
ylabel("Diameter [m]")

F = Isp.*mdot.*g;%Defintion of Isp from lecture7
% It = F*m/mdot
m = It*mdot./F;
V = m.*R.*298./Pc';
figure
plot(Pc,V(:,1))
title("Volume needed to carry required mass of N2 vs Pressure")
xlabel("Tank Pressure [Pa]")
ylabel("Volume [m^3]")

rho =1.021;%g/cm^3
rho=rho/1000;
rho = rho*100^3;
VMMH= m(2)/rho;

disp("--real values--")
disp("Mdot:")
disp(labels'+": "+string(mdot')+"kg/s")
disp("Isp:")
disp(labels'+": "+string(Isp')+"s")
disp("Propellent Mass:")
disp(labels'+": "+string(m')+"kg")
disp("Expansion ratio in both cases is the one found in class: "+string(epsilon))
disp("assuming a density of "+string(rho)+"kg/m^3, the volume of MMH required is: "+string(VMMH)+"m^3")

%% Problem 3
clear
disp("Problem 3")
gam = [1.4,1.25,1.2];
Tc = [300,1300,1650];
R = [297,290,260];
epsilon = 20; %expansion ratio
mdot = .05;%mass flow
labels = ["N2","H2O2","MMH"];
g = 9.8;

Cs = sqrt(gam.*R.*Tc)./(gam.*sqrt((2./(gam+1)).^((gam+1)./(gam-1))));% C star from lecture7
syms Pce
eqn1 = epsilon == (2./(gam+1)).^(1./(gam-1)).*Pce.^(1./gam).*(((gam+1)./(gam-1)).*(1-Pce.^((1-gam)./gam))).^-.5;% definition expansion ratio lecture3
Pce=[vpasolve(eqn1(1)),vpasolve(eqn1(2)),vpasolve(eqn1(3))]; % Pc/Pe
Pce = double(Pce);
Pec = 1./Pce;%Pe/Pc
Cf = sqrt(... % Coeffecient of thrust from lecture7
    ((2*gam.^2)./(gam-1)).*...
    ( (2./(gam+1)).^((gam+1)./(gam-1)) ).*...
    (1-Pec.^((gam-1)./gam)) )+...
    Pec*epsilon;
syms Isp
eqn2 = Cs==Isp*g./Cf;% C star definition from lecture7
Isp = [vpasolve(eqn2(1)),vpasolve(eqn2(2)),vpasolve(eqn2(3))]; 
Isp = double(Isp);
F = Isp.*mdot.*g;%Defintion of Isp from lecture7

disp("Pc/Pe:")
disp(labels'+": "+Pce')
disp("C star:")
disp(labels'+": "+Cs'+"m/s")
disp("Cf:")
disp(labels'+": "+Cf')
disp("Isp:")
disp(labels'+": "+Isp'+"s")
disp("F:")
disp(labels'+": "+F'+"N")
%% 
% PERFORMANCE PARAMETERS N2O2
% 
%  Ae/At                      1.0000   21.805   20.000
%  CSTAR, M/SEC               1542.0   1542.0   1542.0
%  CF                         0.6593   1.7408   1.7297
%  Ivac, M/SEC                1903.0   2828.7   2816.0
%  Isp, M/SEC                 1016.6   2684.3   2667.2
2816.0/g
2667.2/g

 % PERFORMANCE PARAMETERS MMH
 % 
 % Ae/At                      1.0000   24.551   20.000
 % CSTAR, M/SEC               1220.4   1220.4   1220.4
 % CF                         0.6744   1.7695   1.7375
 % Ivac, M/SEC                1514.5   2288.1   2260.6
 % Isp, M/SEC                  823.0   2159.5   2120.4
2260.6/g
2120.4/g

%%
% matlab.codetools.requiredFilesAndProducts("HW4.m")

