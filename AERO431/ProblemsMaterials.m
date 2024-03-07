%% Joshua Oates - A431 - Materials Problem set
clear all
close all
clc

%% stress intensity factor
a = 6; %mm
a = a/1000;% m
% c = 2*a;
sig = 177;%MPa
Y = 1;

K1 = sig*Y*sqrt(pi*a)


% a = 8; %mm
% a = a/1000;% m
% c = 2*a;
% sig = 144;%MPa
% Y = 1;
% 
% K1 = sig*Y*sqrt(pi*a)

%% crack propegations
a = 9; %mm
a = a/1000;% m
% c = 2*a;
sig = 105;%MPa
Y = 1;
crit = 13.6;

K1 = sig*Y*sqrt(pi*a);
K1>crit

%% 6.10

L = 8;%in
A = 12;%in^2
P = 10;%kip
d = .003;%in

E = P*L/(d*A)

%% 6.11
clear all
sigy = 420 %Mpa
sigp = (280+350)/2
epsp = .04+.01*(1/8);
E = sigp/epsp


%% 6.12
clear all
dia1 = 12.7;%mm
L = 50.8;%mm
L = L/1000;
P = 50;%kN
P = P/1000;
dia2 = 12.67494;%mm

dia2 = dia2/1000;
dia1 = dia1/1000;

E = 490/.007;

A = pi*(dia1/2)^2;
sig = P/A;

eps_a = sig/E;
eps_t = (dia1-dia2)/dia1;

nu = eps_t/eps_a

%% 6.13
clear all
dia1 = 12.5;%mm
L = 50;%mm
P = [0
    7
    21
    36
    50
    53
    53
    54
    75
    90
    97
    87.8
    83.3
    ];

d = [0
    .0125
    .0375
    .0625
    .0875
    .125
    .2
    .5
    1
    2.5
    7
    10
    11.5
    ];

A = pi*(dia1/2)^2;
sig = P./A; % kN/mm^2
sig = sig*1000;%MPa
eps = d./L;

sigy = sig(5);
epsy = eps(5);
E = sigy/epsy;

sigr = sig(end);
epsr = eps(end);

[sigu,iu] = max(sig);
epsu = eps(iu);

X = [epsy,epsu,epsr];
Y = [sigy,sigu,sigr];

figure
hold on
plot(eps,sig)
title("Whole Curve")
xlabel("Strain")
ylabel("Stress [MPa]")
scatter(X,Y)


n = 6;
figure
hold on
plot(eps(1:n),sig(1:n))
title("Elastic region")
xlabel("Strain")
ylabel("Stress [MPa]")
scatter(epsy,sigy)

disp("rupture happens at sig = "+string(sigr)+" MPa")
disp("ultimate happens at sig = "+string(sigu)+" MPa")
disp("yeild happens at sig = "+string(sigy)+" MPa")
disp("youngs modulous is: "+string(E)+" MPa")

