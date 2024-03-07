%% Exam 1 - A355 - question 5 - Joshua Oates
clear all
close all
clc

g = 9.81; % m/s^2
R = 8.314; % J/(mol*K)
avo = 6.0221408e+23; %atoms / mol
k = 1.38e-23; %J/K

P0 = 1; %Pa
h0 = 80; % km
T0 = -80;% C
T0 = T0 + 273.15; %K
L = 3; % K/km
M = 16; % g/mol
M = M/1000; % kg/mol



ng0 = P0/(k*T0); %atoms/m^3
molm = ng0/avo; %mol/m^3

rho0 = molm*M; %kg/m^3

alt = linspace(80,800)'; % km
T = T0+L*(alt-h0); % K

H = (R*T)./(M*g); % m
H = H/1000; % km

rho = rho0*exp(-(alt-h0)./H);


%%%%%%%%%%%%%%%%%%%%%
figure
hold on
plot(alt,rho)
gca().set('YScale','log')
xlabel('altitude [km]')
ylabel('particle density [kg/m^3]')
title('Particle Density Vs Altitude')


% part b

m = 1; %kg
A1 = 100; %cm^2
A1 = A1/(100^2); % m^2
A2 = 10; %m^2
A = linspace(A1,A2)';

z = 800; % km
Cd = 2.2;

r_e = 6378;%km
mu_e = 3.98600e14;%m^3/s^2
a = r_e+z; % km
a = a*1000; % m

da_dt = -rho(end).*(Cd.*A./m).*(mu_e*a)^.5;% m/s
da_dt = da_dt/1000; %km/s

% da_dt2 = -rho*(Cd*A2/m)*(mu_e*a)^.5;% m/s
% da_dt2 = da_dt2/1000; %km/s

Hmax = H(end);


L = -Hmax./da_dt;
L = L./(60*60*24); % days
L1 = L(1)

%%%%%%%%%%%%%%%%%%%%%%
figure 
hold on 
plot(A,L)
gca().set('YScale','log')
gca().set('XScale','log')
xlabel('area [m^2]')
ylabel('Snapshot Lifetime [days]')
title('Lifetime Vs Area')


% disp
disp("The average mass density at 800 km is: "+string(rho(end))+" kg/m^3")
disp("The snapshot lifetime of the satalite at 800 km with a cross sectional area of 100 cm^2 is: "+string(L(1))+" days.")
disp("The snapshot lifetime of the satalite at 800 km with a cross sectional area of 10 m^2 is: "+string(L(end))+" days.")
disp("This seems impossibly short to me, perhaps in my calculations I had a units error or perhaps the numbers I used as input are not quite realistic. I suspect that somewhere there is a units error but after significant trouble shooting I couldn't find it. At any rate, It can be seen that as expected the larger area there will have a shorter lifespan than the smaller one.")






