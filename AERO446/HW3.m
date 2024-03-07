clear all
close all
clc

%% P1

mu_e = 398600;%km^3/s^2
syms r
T = 6*3600;
eqn = T == 2*pi*r.^(3/2)./sqrt(mu_e); % s
r = double(solve(eqn,r));


% r = 8000;%km
r_e = 6378;%km

T = 2*pi*r.^(3/2)./sqrt(mu_e); % s

% r_e = r*sin(theta_shadow)
theta_shadow = asin(r_e/r);
theta_shadow = theta_shadow*2;

percent_eclipse = theta_shadow/(2*pi);

disp("P1")
T_eclipse = T*percent_eclipse
%%%%%%confirmed

%% P2
mu_e = 398600;%km^3/s^2
r = 8000;%km
r_e = 6378;%km

T = 2*pi*r.^(3/2)./sqrt(mu_e); % s

% r_e = r*sin(theta_shadow)
theta_shadow = asin(r_e/r);
theta_shadow = theta_shadow*2;

percent_eclipse = theta_shadow/(2*pi);

disp("P2")
T_eclipse = T*percent_eclipse
%%%%%%confirmed


%% P3
P_req = 500;% W
dod = .5; % depth of discharge
eta =1; % no reason not to assume good efficiency

disp("P3")
E_req = P_req*T_eclipse/(dod*eta); % J
E_req = E_req*0.000277778 % wh



%% P4
P_A = 245; %W/m^2
syms A
P_generated = P_A*A*(1-percent_eclipse);

disp("P4")
A = double(solve(P_generated == P_req))

%% P5
P_e = P_req+100;
P_ave = percent_eclipse*P_e+(1-percent_eclipse)*P_req;
syms A
P_generated = P_A*A*(1-percent_eclipse);

disp("P5")
A = double(solve(P_generated == P_ave))

%% P6
% P_e = P_req+100;
P_ave = percent_eclipse*P_e+(1-percent_eclipse)*P_req;
syms A t

E_generated = (int(sin(t),[0 pi])/pi)*A*P_A;
P_generated = E_generated/2;

disp("P6")
A = vpa(solve(P_generated == P_ave))

%% P7

disp("P7")
effscale10 = 1*(1-.02)^10;
P_ave = percent_eclipse*P_e+(1-percent_eclipse)*P_req;
syms A t

E_generated = (int(sin(t),[0 pi])/pi)*A*P_A*effscale10;
P_generated = E_generated/2;

A = double(solve(P_generated == P_ave))


%% P8

disp("P8")
% P1 repeat
syms r
T = 6*3600;
eqn = T == 2*pi*r.^(3/2)./sqrt(mu_e); % s
r = double(solve(eqn,r));

mu_e = 398600;%km^3/s^2
% r = 8000;%km
r_e = 6378;%km

T = 2*pi*r.^(3/2)./sqrt(mu_e); % s

% r_e = r*sin(theta_shadow)
theta_shadow = asin(r_e/r);
theta_shadow = theta_shadow*2;

percent_eclipse = theta_shadow/(2*pi);
T_eclipse = T*percent_eclipse;

%%%%%%%%%%

effscale10 = 1*(1-.02)^10;
P_ave = percent_eclipse*P_e+(1-percent_eclipse)*P_req;
syms A t

E_generated = (int(sin(t),[0 pi])/pi)*A*P_A*effscale10;
P_generated = E_generated/2;

A = vpa(solve(P_generated == P_ave))


