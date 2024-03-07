%% Joshua Oates final
close all
clear all
clc

%% Part 3
mu = 398600;
Re = 6378;
kb = 1.38e-23;
m_e = 9.11e-31; %kg
q = 1.602e-19; %C

radius = .5; 
z = 300;

Te = 1500;
eta = 5e11;

A = pi*radius^2;

v_sc = sqrt(mu/(Re + z))*1e3; %m/s

I_e = @(V) .25*q*eta*sqrt((8*kb*Te) / (pi*m_e)) * exp(q*V/(kb*Te));
I_i = q*eta*v_sc*A;

Vf = fzero(@(V) I_e(V)-I_i, 0)
I_e = I_e(linspace(-1,0));

figure
hold on 
plot(I_e)
yline(I_i)

xlabel("S/C Voltage")
ylabel("I [A]")
legend('electron', 'ion')

mu_0 = 4*pi*1e-7;

f_p = 8.98*sqrt(eta)
l_db = sqrt(mu_0*kb*Te / (eta*q^2))
lam = (4*pi*eta*l_db^3)/3

set(gca,"XScale","log")
set(gca,"YScale","log")

%% Part 4
eta = 95e6;
M = 1.67e-27;
k = 2;
B = 30500e-9; % Tesla
mu = pi*4e-7;
radius = 4.3;
v = sqrt( (k^2*B^2) / (radius^6*eta*mu*M) ); % m/s
v = v*1e-3; % km/s

disp("The velocity of the solar wind is " + v + " km/s")

%% Part 5
clear
d = readmatrix("data.csv");
z = d(:,1);
dz = mean(diff(z));
dz = dz*1000;
ne = d(:,2);
TEC = cumsum(ne.*dz);

fsig = linspace(1e6,10e9);

k2 = 40.31;
xr = (k2./(fsig.^2)).*TEC;
[X,Y] = meshgrid(fsig,z);
figure
hold on
surf(X,Y,xr)
xlabel("Frequency")
ylabel("z")
zlabel("Extra Range")

set(gca,"XScale","log")
% set(gca,"YScale","log")
set(gca,"ZScale","log")

figure
hold on
alts = [200 400 600 1000 2000];
for i = alts
    xr = (k2./(fsig.^2)).*TEC(z==i);
    plot(fsig,xr)
end

set(gca,"XScale","log")
set(gca,"YScale","log")
legend(string(alts)+" km")
xlabel("Frequency")
ylabel("Extra Range")



%% Part 6
d = 1.25; % m
T0 = 55; % K temp O2
sig = 5.6704e-8; % steffan boltzmann
SolarConst = 1.361e3;% W/m^2
TDS = 2.75;% K temp deep space

As = 4*pi*(d/2)^2; % area of sphere
Ac = pi*(d/2)^2; % cross sectional area of sphere

% epsilon = alpha

% e_m = .76; %mylar
% e_k = .72; %kapton
% ek = e_m; % interior layer
% ei = e_k; % inner layer
% eo = e_m; % outer layer
% %   layers:
% %          k | m t m t m h m t m t | m
% %          i | n   n   n   n   n   | o
%
% n = 5;
% e_eff = (ei^-1+eo^-1-1+2*n*ek^-1-n)^-1 % effective

eps_o = .1;  % outer
eps_m = .07; % interior (middle)
eps_i = .2;  % inner

n = 20;
eps_eff = (eps_i^-1 + eps_o^-1 -1 + 2* n* eps_m^-1 -n)^-1; % effective

Qin = SolarConst*Ac*eps_eff;
Qout = sig*As*eps_eff*(T0^4 - TDS^4);
Qnet = Qin-Qout;
disp(n+" interior layers will give a net heat transfer of "+Qnet+" watts being transfered to the LOX")

syms Tos % Toutside
eqn = Qnet == eps_o*sig*As*(Tos^4 - TDS^4);
Tos = solve(eqn,Tos);
Tos = double(max(real(Tos)));
disp("The average outter layer temperature must be "+Tos+" K to have enough exitance leaving the S/C to maintain 3W total gained heat. This calc is based on the emittance of the outer most layer of the MLI exclusively")

