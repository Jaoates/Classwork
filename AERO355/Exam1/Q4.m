%% Exam 1 - A355 - question 4 - Joshua Oates
clear all
close all
clc
d = load("AtmosData.txt");

%% constants
k = 1.38e-23; %J/K
g = 9.81; % m/s^2
R = 8.314; % J/(mol*K)
avo = 6.0221408e+23; %atoms / mol

% reference pressure
P0 = 1; % at 80 km
h0 = 80; % km
T0 = 180; %k


%%
n = length(d);

% create temperature vectors for solar min and max
Tmin = linspace(T0,700,n);
Tmax = linspace(T0,1800,n); % K

% breifely to vectors
alt = d(:,1);
Mmin = d(:,2);
Mmax = d(:,3);
Amin = d(:,4);
Amax = d(:,5);

Aave = (Amin+Amax)/2;
Tave = (Tmin+Tmax)/2;
Mave = (Mmin+Mmax)/2;

clear d
d = struct('alt',alt,'Mmin',Mmin,'Mmax',Mmax,'Mave',Mave,'Amin',Amin,'Amax',Amax,'Aave',Aave,'Tmin',Tmin','Tmax',Tmax','Tave',Tave');
clear alt Mmin Mmax Amin Amax Tmin Tmax Tave Mave Aave

%% calcs
d.LapseMax = (max(d.Tmax)-T0) / (max(d.alt)-h0); % K/km
d.LapseMin = (max(d.Tmin)-T0) / (max(d.alt)-h0);
d.LapseAve = (max(d.Tave)-T0) / (max(d.alt)-h0);

d.Hmax = (R*d.Tmax)./(d.Mmax*g); % m
d.Hmax = d.Hmax/1000; % km
d.Pmax = P0*exp(-(d.alt-h0)./d.Hmax);

d.Hmin = (R*d.Tmin)./(d.Mmin*g); % m
d.Hmin = d.Hmin/1000; % km
d.Pmin = P0*exp(-(d.alt-h0)./d.Hmin);

d.Have = (R*d.Tave)./(d.Mave*g); % m
d.Have = d.Have/1000; % km
d.Pave = P0*exp(-(d.alt-h0)./d.Have);


d.ngmax = d.Pmax./(k*d.Tmax);
d.ngmin = d.Pmin./(k*d.Tmin);
d.ngave = d.Pave./(k*d.Tave);

d.massMax = d.Mmax/avo; % kg
d.massMin = d.Mmin/avo;
d.massAve = d.Mave/avo;

d.veloMax = sqrt(8*k*d.Tmax./(pi*d.massMax)); % m/s
d.veloMin = sqrt(8*k*d.Tmin./(pi*d.massMin));
d.veloAve = sqrt(8*k*d.Tave./(pi*d.massAve));

d.sigMax = pi*(d.Amax);
d.sigMin = pi*(d.Amin);
d.sigAve = pi*(d.Aave);

d.mfpMax = 1./(4*d.ngmax.*d.sigMax); 
d.mfpMin = 1./(4*d.ngmin.*d.sigMin); 
d.mfpAve = 1./(4*d.ngave.*d.sigAve); 


% disp(d.Pmin)
% d.Pmin = P0*(1+(d.LapseMin./T0)*(d.alt-h0)).^-(d.Mmin*g/(1000*R*d.LapseMin));
% disp(d.Pmin)


%% plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pressure
figure
hold on
plot(d.alt,d.Pmin)
plot(d.alt,d.Pmax)
plot(d.alt,d.Pave)
xlabel('altitude [km]')
ylabel('pressure [Pa]')
title('Pressure Vs Altitude')

% gca().set('XScale','log')
gca().set('YScale','log')

legend('Solar Min','Solar Max','Solar Average','location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% particle density
figure
hold on
plot(d.alt,d.ngmin)
plot(d.alt,d.ngmax)
plot(d.alt,d.ngave)
xlabel('altitude [km]')
ylabel('density [1/m^3]')
title('Particle Density Vs Altitude')

% gca().set('XScale','log')
gca().set('YScale','log')

legend('Solar Min','Solar Max','Solar Average','location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermal Velo
figure
hold on
plot(d.alt,d.veloMin)
plot(d.alt,d.veloMax)
plot(d.alt,d.veloAve)
xlabel('altitude [km]')
ylabel('thermal velocity [m/s]')
title('Thermal Velocity Vs Altitude')

% gca().set('XScale','log')
% gca().set('YScale','log')

legend('Solar Min','Solar Max','Solar Average','location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MFP
figure
hold on
plot(d.alt,d.mfpMin)
plot(d.alt,d.mfpMax)
plot(d.alt,d.mfpAve)
xlabel('altitude [km]')
ylabel('mean free path [m]')
title('Mean Free Path Vs Altitude')

% gca().set('XScale','log')
gca().set('YScale','log')

legend('Solar Min','Solar Max','Solar Average','location','best')

%% discussion
disp("The following data is as calcuated at 400 km:")
disp("  Pressure in Pa for solar min max and average in that order: ")
disp(string([d.Pmin(end);d.Pmax(end);d.Pave(end)]))
disp("  Particle density in atoms/m^3 for solar min max and average in that order: ")
disp(string([d.ngmin(end);d.ngmax(end);d.ngave(end)]))
disp("  Thermal velocity in m/s for solar min max and average in that order: ")
disp(string([d.veloMin(end);d.veloMax(end);d.veloAve(end)]))
disp("  Mean free path in m for solar min max and average in that order: ")
disp(string([d.mfpMin(end);d.mfpMin(end);d.mfpMin(end)]))

disp("From the plots you can see basic trends in the data. As expected particle density and pressure both steeply decrease with altitude and in both cases they decrease the least at solar max. Thermal Velocity and mean free path as expected increase with altitude. This makes intuitive sense because if there are less particles there will be more room between them to move before experiencing a collision.")

%% units check
% clc
% 
% u = symunit;
% g = u.m/u.s^2;
% R = u.J/(u.mol*u.K);
% alt = u.km;
% M = u.kg/u.mol;
% T = u.K;
% L = u.K/u.km;
% 
% P = u.Pa;
% k = u.J/u.K;
% 
% ng = P/(k*T);
% disp(simplify(ng))


% x = M*g/(1000*R*L);
% disp(simplify(x))
% 
% H = R*T/(M*g);
% disp(simplify(H))
% H = 1000*R*T/(M*g);
% disp(simplify(H))


