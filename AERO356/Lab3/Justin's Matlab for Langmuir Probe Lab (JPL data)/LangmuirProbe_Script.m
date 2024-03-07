%{
Justin Self
Cal Poly
Space Environments II Lab
May 22 2023
LAB 3: LANGMUIR PROBE
%}

clear all; close all; clc;

%% Knowns and Extracting data

% Knowns
% Langmuir probe dimensions
d = 0.02; % inches
l = 0.075; % inches
area_inches = pi*(d/2)^2 + l*(pi*d);
insq2m2 = 0.00064516; % in2 to m2
Aprobe = insq2m2 * area_inches; % in m2

disp("-------------- JPL DATA --------------")
JPL.raw = readmatrix("LangmuirProfessional.xlsx");
JPL.v = JPL.raw(:,1);                                   % v
JPL.i = JPL.raw(:,2);                                   % amps
x = JPL.v;
y = JPL.i;

% FOR STUDENT DATA, UN-COMMENT THIS PORTION BELOW:

%{
% Student Data
disp("-------------- STUDENT DATA --------------")
student.raw = readmatrix("myfilename.xlsx");
student.v = student.raw(:,1);                                   
student.i = student.raw(:,2);                                   
x = student.v;
y = student.i;
%}

% Plot raw data
figure()
plot(x,y,'.')
xline(0)
yline(0)

% Graph Raw Data 
ylim padded 
xlim tight 
xLab = xlabel('Voltage [V]','Interpreter','latex'); 
yLab = ylabel('Current [A]','Interpreter','latex'); 
plotTitle = title('JPL Langmuir Probe (Raw Data)','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 

% Shift plot up so non-negative
minVal.x = abs(min(x));
minVal.y = abs(min(y));

% Extract roughly linear portion by eye
xlin = find(x>0 & x<30);
ylin = y(xlin);                                         % <- linear y vector (rough)
xlin = x(xlin);                                         % <- linear x vector (rough)


% Create shifted vectors (nonnegative)
shifted.x = x + minVal.x;
shifted.y = y + minVal.y;
shifted.xlin = xlin + minVal.x;
shifted.ylin = ylin + minVal.y;

% Plot raw shifted
figure()
semilogy(shifted.x,shifted.y,'.')
hold on
semilogy(shifted.xlin,shifted.ylin,'r','linewidth',2)

ylim padded 
xlim tight 
xLab = xlabel('Voltage [V]','Interpreter','latex'); 
yLab = ylabel('Current [A]','Interpreter','latex'); 
plotTitle = title('JPL Langmuir Probe (Shifted) ','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Raw data','Rough Linear Portion','interpreter','latex','Location', 'best')

%% Analysis

%{
Things we want to calculate: 
Te = electron temperature
np, ni, ne = density of plasma, ions, electrons
phi_p = plasma potential
Ii = saturation current ion
Ie = saturation current electron
debye_length 
wp = plasma frequency
Vf = floating potential (V when Ii = Ie)
Vi,th = ion thermal velocity

%}

% Knowns
k = 1.380649e-23;                           % Boltzman constant; m2 kg s-2 K-1

% Ar atomic mass = 39.95
amus = 1.660538921e-27;                     % atomic mass units in a kg
argon_mass = 39.95 * amus;                  % mass of an argon atom; kg
qe = 1.60217663e-19;                        % charge of an electron, C
me = 9.1093837e-31;                         % mass of an electron, kg

% Find electron temperature, Te, in eV
deltaV = shifted.xlin(1) - shifted.xlin(end);
ratioI = shifted.ylin(1) ./ shifted.ylin(end);

Te_k = qe*deltaV / (k*log(ratioI)); % in K

Te_eV = deltaV / log(ratioI);    % in eV now      

disp("Floating potential (Vf at 0 amps) is: 0 ")

disp("Electron temperature is: " + Te_eV + " eV")
disp("Electron temperature is: " + Te_k + " K")

% Find Ve,th == electron thermal velocity
Ve_th = sqrt(8*k*Te_k / (pi*me) );

disp("Electron thermal velocity is: " + Ve_th/1000 + " km/s")

%% Saturation currents

% Look at original plot; find where left side of plot dI/dV = 0.
x_for_Iis = find(shifted.x > 10 & shifted.x < 15);
y_for_Iis = shifted.y(x_for_Iis);
Iis = mean(y_for_Iis);

%Iis = shifted.y(1);
disp("Ion (Ar) saturation current is: " + Iis + " amps")

% Ion density
ni = (Iis/(qe*0.6)) * (k * Te_k / argon_mass)^(-.5) * (1/Aprobe);
disp("Ion density (Ar) is: " + round(ni,4) + " particles/m3")

% Electron saturation current
x_for_Ies = find(shifted.x > 90 & shifted.x < 100);
y_for_Ies = shifted.y(x_for_Ies);
Ies = mean(y_for_Ies);

% Ies/Iis ~180 CHECK
disp(" ")
disp("      CHECK: Ies/Iis is: " + Ies/Iis)
disp("      Close-ish to 180!")
disp(" ")

disp("Electron saturation current is: " + Ies + " amps")

% Electron density
ne = 4  * Ies / (qe * Ve_th * Aprobe);
disp("Electron density is: " + round(ne,4) + " particles/m3")

% Check ne ~ ni
disp(" ")
disp("      CHECK: ne/ne (should be close to 1) is: " + ne/ni)
disp("      Close-ish to 1!")

% Debye length (optional)
epsilon0 = 8.85e-12; % permittivity of free space
lambda = sqrt(  (epsilon0 * k * Te_k) / (ne * qe^2)   );

disp(" ")
disp("Debye Length is: " + lambda + " m")

%% Plot shifted with analysis lines

% Need linear portion fit line
coeffs = polyfit(shifted.xlin,log(shifted.ylin),1);
p = shifted.x .*coeffs(1) + coeffs(2);

% top (Ies) intersection 
int.Ies = [66.174965,Ies];

% intersection of int.Ies(x) and Iis (horizontal line)
plasmaPotential = [66.174965,4.6420455e-4];

% Plot it
figure()
semilogy(shifted.x,shifted.y,'.')
% plot(shifted.x,shifted.y,'.')
hold on
yline(Ies,'r--')
yline(Iis,'b--')
semilogy(shifted.x,exp(p),'r','linewidth',1)
% plot(shifted.x(1:600),exp(p(1:600)),'r','linewidth',1)

p2 = plot(plasmaPotential(1),plasmaPotential(2),'r','linewidth',2);
p2.Marker = "diamond";

% plot vertical line
xline(int.Ies(1))

ylim padded 
xlim tight 
xLab = xlabel('Voltage [V]','Interpreter','latex'); 
yLab = ylabel('Current [A]','Interpreter','latex'); 
plotTitle = title('JPL Langmuir Probe (Analysis) ','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Raw data','$I_{es}$','$I_{is}$','Slope of linear portion (log plot)','$\phi_p$','interpreter','latex','Location', 'best')

 %% Discussion
 disp(" ")
 disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
 disp("Discussion")
 disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
 disp(" ")
 disp("Q:   What are the plasma parameter values?")
 disp("         A: See above outputs ")
 disp("Q:   Why are these plasma parameters important?")
 disp("         A: ______ ")
 disp("Q:   How does the student data compare to the professional data?")
 disp("         A: ______ ")
 disp("Q:   What are the sources of error?")
 disp("         A: (methods of extracting data; log plots; arbitrary choices on linear portion; subjective.) ")