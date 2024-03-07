%% 446 - HW1 - Joshua Oates
clear all
close all
clc
addpath("C:\joshFunctionsMatlab\")

%% orbital period and FOR
clear

SC = ["WV3";"1W";"O3b";"GEO"];

mu_e = 398600;%km^3/s^2
r_e = 6378;%km

Alt = [650;1070;8000;35786];
r = r_e+Alt;

T  = 2*pi*r.^(3/2)./sqrt(mu_e); % s
T = T/60; %min

FOR = 2*asind(r_e./(r));

table(SC,T,FOR)


%% Columbia
clear all

C_c = 30200;
r_c = C_c/(2*pi);
r_e = 6378;
V_e = (4/3)*pi*r_e^3;
V_c = (4/3)*pi*r_c^3;

Vratio = V_c/V_e;

mu_e = 398600;%km^3/s^2
mu_c = mu_e*Vratio;

syms r
eqn = 24*60*60 == 2*pi*r.^(3/2)./sqrt(mu_c); % s
r = double(solve(eqn,r));
z = r-r_c;
FOR = 2*asind(r_c./(r))

%% WV3
% FOV
clear all
mu_e = 398600;%km^3/s^2

z = 650;
r_e = 6378;
r = r_e+z;

FOV = 5;
FOV = 2*atand(FOV/z)

% view time
x = 10;
v = sqrt(mu_e/r);
T = x/v



