
%% 0 
close all;
clear all;
clc

addpath("C:\joshFunctionsMatlab\")

%% toolbox
clear all;
DM = importdata("DMout1.mat");
R = 287;% J/kg.K

n = 100;
M = linspace(0,1,n);

Tr = zeros(1,n);
Pr = zeros(1,n);
rhor = zeros(1,n);

for i = 1:n
    [Tr(i),Pr(i),rhor(i)] = joshIsentropicToolbox(M(i));
end

figure
hold
plot(M,Tr.^-1)
plot(M,Pr.^-1)
plot(M,rhor.^-1)
legend("Tr","Pr","rhor")

gamma =1.4 ; 
a =@(T) sqrt(gamma*R*T);
Ma = @(U,T) U/a(T);
Ttotal = DM(1).T_amb;
Ptotal = DM(1).Amb_Press;
rhototal = 1.204; % kg/m^3
Cp = 1005 ;%J/kg.k

for i=1:length(DM)
    for j = 1:length(DM(i).V_fromA)
        DM(i).Ma(j) = Ma(DM(i).V_fromA(j),Ttotal);
    end
end

for i = 1:length(DM)
    for j = 1:length(DM(i).Ma)
        [DM(i).Tr(j),DM(i).Pr(j),DM(i).rhor(j)] = joshIsentropicToolbox(DM(i).Ma(j));
        DM(i).T_isen(j)=Ttotal/DM(i).Tr(j);
        DM(i).P_isen(j)=Ptotal/DM(i).Pr(j);
        DM(i).rho_isen(j)=rhototal/DM(i).rhor(j);
        DM(i).error(j) = (DM(i).P(j)-DM(i).P_isen(j))/DM(i).P_isen(j);
        DM(i).ds(j) = Cp*log(DM(i).Tr(j))-R*log(DM(i).Pr(j));
    end
end

disp("------------ISF------------")
disp("My work has the following results:")
disp("I calculated that as compared to the measured values of P, isentropic relations find the following errors as a percent and delta s in kJ/kg.K in that order at each station:")
disp("run1:")
disp("error:")
disp(DM(1).error*100)
disp("delta s:")
disp(DM(1).ds)

disp("run2:")
disp("error:")
disp(DM(2).error*100)
disp("delta s:")
disp(DM(2).ds)

disp("run3:")
disp("error:")
disp(DM(3).error*100)
disp("delta s:")
disp(DM(3).ds)

disp("run4:")
disp("error:")
disp(DM(4).error*100)
disp("delta s:")
disp(DM(4).ds)

disp("This loss of entropy generally less than a percent seems to be plenty within limits to call the windtunnel isentropic. Most likely the flow in the center of the tunnel is closer to isentropic flow than that on the edges, since we are assumming that the loss of entropy comes mostly from friction with the walls and possibly some leaking in the wind tunnel which will cause a loss of total pressure resulting in a measured loss of entropy.")

disp("I calculated the following density ratios based on Isentropic assumptions:")
disp("run1:")
disp("rhoTotal/rho:")
disp(DM(1).rhor)

disp("run2:")
disp("rhoTotal/rho:")
disp(DM(2).rhor)

disp("run3:")
disp("rhoTotal/rho:")
disp(DM(3).rhor)

disp("run4:")
disp("rhoTotal/rho:")
disp(DM(4).rhor)

disp("These values seem to be close to the plots of rho ratio vs Ma, the highest Ma is around .1 and by the graph this should yeild a rho ratio of about 1/.99 or 1.005. With such low mach numbers it is slightly hard to tell however, rho almost doesn't change and incompressible is probably not a bad assumption for this process in general.")
%% FK
% FK was completed by hand, see PDF

%% CM1
clear all
Ue = 15;%m/s
L = 10;%m
Xdot = @(t,X) X*Ue/L;

dbulb = .15;
dcond = .1;
V0 = (4/3)*pi*(dbulb/2)^3;
P = 199325;
rho = .9108;
m = rho*V0;

tspan = [0 .75];

options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
X0 = 0;
[t1,X1] = ode45(Xdot,tspan,X0,options);
X0 = 3;
[t2,X2] = ode45(Xdot,tspan,X0,options);


rho = @(v) m/(V0+v);
v = @(x) pi*(dcond/2)^2*x;

for i = 1:length(X2)
    rhov(i) = rho(v(X2(i)));
end

figure
plot(t2,rhov)
title("rho vs time")
xlabel("t [s]")
ylabel("density [kg/m^3]")





%% CM2
clear all

syms U0 d rho% U0 is [m/s] and d is [m] but i will consider they themselves to be units of length and velocity respectively
A0 = pi*(d/2)^2; % [d^2]
vdot0 = U0*A0; % mass flow rate at exit in terms of U0 and d
mdot0 = vdot0*rho;

% in this case below, x means times or function handle name
u_xdxU0=@(x,r)  (5/x)*exp(-50*r^2/x^2); % a function that returns U/(d*U0) at any x and r (used for when U0 and d are not known)
% true u(x,r) given by (5/x)*exp(-50*r^2/x^2)*U0*d

n = 100;
x = 10; % units of [d]
r = linspace(0,5,n); % units of [d]
Uprofile = zeros(1,n);

for i = 1:n
    Uprofile(i) = u_xdxU0(x,r(i)); % units of [(m/s)/(d*U0)] ie non dimensional
end

figure
plot(r,Uprofile)
xlabel("r [d]")
ylabel("U [U0]")

clear Uprofile r x n


vdot_xdxU0 =@(x) (1/10)*pi*x; % gives vdot in units of [d^3/s] takes x [d] 
% true vdot given by vdot = (1/10)*d*pi*U0*x

mdot_xdxU0 =@(x) vdot_xdxU0(x)*rho;

n = 100;
x = linspace(0,5,n); % units of [d]
vdotprofile = zeros(1,n);

for i = 1:n
    vdotprofile(i) = vdot_xdxU0(x(i)); % units of [d^3] 
end
figure
plot(x,vdotprofile)
xlabel("x [d]")
ylabel("vdot [d^2.U0]")

% volumetric flow rate is effectively a velocity times an area, to get
% average velocity over an area we need just divide by that area. Up until
% now we have been integrating by infinite bounds but we now must choose a
% control area lets assume that area in question is the circlular cross
% section of the cone with angle 11.8 degrees at the x we care about
theta = 11.8; % degrees
r = @(x) sind(theta)*x; % [d]
A = @(x) pi*r(x)^2; % [d^2]
aveVelo_xdxU0 = @(x) vdot_xdxU0(x)/A(x); % [U0]
% true aveVelo given by U0*d*vdot_xdxU0(x)/A(x)

aveVelo_xdxU0_10=aveVelo_xdxU0(10);

n = 100;
x = linspace(5,12,n); % units of [d]
aveVeloprofile = zeros(1,n);
for i = 1:n
    aveVeloprofile(i) = aveVelo_xdxU0(x(i)); % units of [d^3] 
end
figure
plot(x,aveVeloprofile)
xlabel("x [d]")
ylabel("Average Velocity [U0]")

clear d rho U0 A0
d = 5;%cm
d = d/100;%m
U0 = 100;%m/s
A0 = (pi/4)*d^2;
[~,~,rho] = joshStdAtm(); % assume rho is at sea level
mdot = @(x) vdot_xdxU0(x)*rho*U0*d;
mdot10d = mdot(10*d);
mdot0 = U0*A0*rho;

% conservation of mass is not broken, the additional mass in the flow is
% accounted for by air which is accelerated by the initial flow.

disp("------------CM2------------")
disp("My work has the following results:")
disp("Mass flow rate [kg/s] at the nozzle is: "+string(mdot0))
disp("Average Velocity [m/s] of flow at 10d is: "+string(aveVelo_xdxU0_10))
disp("Mass flow rate [kg/s] at 10d is: "+string(mdot10d))
disp("there is a linear increase in mass flow rate as distance increases from the nozzle, but this is to be expected as jet entrainment will cause nearby stationary fluid to become accelerated and enter the jet. this exchange of momentum causes the velocity of the jet to drop but increases the area. This is only possible in a continuum with a sufficeintly viscous fluid.")
disp("I chose to only plot one half of the velocity profile, but it is mirrored across the U axis to create a symetric profile.")
disp("The integrating bounds for all integrals are infinity becuase there is not a sharp cuttoff to jet entrainment, rather it decreases radially from the jet and asyptotically approaches zero. This means however that theoretically at least, any arbitrarily far away particle will have some small change in momentum from this jet. For that reason it is logical to integrate over indefininte bounds.")
disp("I found vdot by inegrating U(r) over an infinite cylindrical area and then found average velocity U-bar by dividing by the crossectional area of the cone at the x point of interest.")