%% depth of errosion

t = 3;%y
dm_dt = 38.3;%g/day
rho = 7.3;%g/cm^3
A = .4;%m^2

y = 365.25;% days/y

rho = rho*(1e6);% rho in meters
t= t*y;% days
dm = dm_dt*t; % total mass loss

depth = dm/(A*rho);% in meters
depth = depth*1000 % in mm

%% orbital decay expected life
da_dt = -5.3; %mm/s
da_dt = da_dt/1000; %m/s
T = 256;%K
M = 18.2; %g
M = M/1000; %kg
R = 8.314;
g = 9.8; %m/s^2

H = R*T/(M*g);

life = -H/da_dt;
life = life/(60*60*24)

%
% da_dt = -8.3; %mm/s
% da_dt = da_dt/1000; %m/s
% T = 431;%K
% M = 12.6;
% M = M/1000
% R = 8.314;
% g = 9.8; %m/s^2
% 
% H = R*T/(M*g);
% 
% life = -H/da_dt;
% life = life/(60*60*24)


%% obital decay rate
clear

A = 8; % m^2
m = 214; % kg
z = 322; % km
Cd = 2.2;
rho = 1e-11; % kg/m^3

mu_e = 3.98600e14;%m^3/s^2

r_e = 6378;%km

a = r_e+z; % km
a = a*1000; % m

da_dt = -rho*(Cd*A/m)*(mu_e*a)^.5;% m/s

da_dt = da_dt*1000 %mm/s

%% Brightness
h = 332;
B = 10^(7-0.0129*h)

%% Sputtering
A = 3.1;% m^2
fluxi = 9.1e22; 
mt = 6.7e-26;%kg
mt = mt*1000;%g

t=1.9;%y
y = .48;

fluxsp = fluxi*y;

dm = fluxsp*A*t*mt

% %%
% A = 6.6;% m^2
% fluxi = 7.9e22; 
% mt = 6.7e-26;%kg
% mt = mt*1000;%g
% 
% t=1.8;%y
% y = .21;
% 
% fluxsp = fluxi*y;
% 
% dm = fluxsp*A*t*mt

%% mass loss AO
t = 1.9;%y
A = 3.1;%m^2
A = A*10000;%cm^2
rho = 8.1;%g/cm^3


flux = 9.5e23;
flux = flux/(10000);
E =4.9e-24;

dm = rho*E*flux*A*t

