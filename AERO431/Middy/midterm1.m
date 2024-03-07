%% Joshua Oates - Midterm 1

%% Q11 6.3
a = 3; %mm
a = a/1000;% m
% c = 2*a;
sig = 148;%MPa
Y = 1;

K1 = sig*Y*sqrt(pi*a)

%% Q20 2.1

clear all
close all
clc

w1 = 29;%mm
w1 = w1/1000;%m
w2 = 79;
w2 = w2/1000;

Abar = w2*w1;
t = 9;%mm
t=t/1000;
tau_allow = 88;%MPa
tau_allow = tau_allow*1e6;

T_allow = 2*Abar*t*tau_allow

L = 4;
G = 75e9;%Pa

phi = (T_allow*L/(4*Abar^2*G*t))*(2*(w1+w2))
rad2deg(phi)

%% Q21 2.1b

clear all
close all
clc

w1 = 27;%mm
w1 = w1/1000;%m
w2 = 78;
w2 = w2/1000;

Abar = w2*w1;
t = 10;%mm
t=t/1000;
tau_allow = 88;%MPa
tau_allow = tau_allow*1e6;

T_allow = 2*Abar*t*tau_allow

L = 4;
G = 75e9;%Pa

Tactual = 3763

phi = (Tactual*L/(4*Abar^2*G*t))*(2*(w1+w2))
rad2deg(phi)

%% Q22 2.2
clear all
close all
clc

T = 55;

w1 = 28;%mm
w1 = w1/1000;%m
w2 = 74;
w2 = w2/1000;

Abar = w2*w1;
t = 8;%mm
t=t/1000;

tau_ave = T/(2*t*Abar)

%% Q23 2.6
clear all

b = 5;%mm
A = 188;
C = 76;
B = 310;

b = b/1000;%m
A = A/1000;
C = C/1000;
B = B/1000;



h = [2*A,2*C];

I = (1/12)*b.*h.^3;

syms x

x = vpa(solve(I(1)*x == I(2)*(B-x)))*1000



