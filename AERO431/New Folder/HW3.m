%% Problems with thin walled sections

%% 2.1
clear all
close all
clc

w1 = 30;%mm
w1 = w1/1000;%m
w2 = 70;
w2 = w2/1000;

Abar = w2*w1;
t = 10;%mm
t=t/1000;
tau_allow = 80;%MPa
tau_allow = tau_allow*1e6;

T_allow = 2*Abar*t*tau_allow

L = 4;
G = 75e9;%Pa

phi = (T_allow*L/(4*Abar^2*G*t))*(2*(w1+w2))
rad2deg(phi)

%% 2.2
clear all
close all
clc

T = 50

w1 = 30;%mm
w1 = w1/1000;%m
w2 = 70;
w2 = w2/1000;

Abar = w2*w1;
t = 10;%mm
t=t/1000;

tau_ave = T/(2*t*Abar)

%% 2.3
I = (1/12)*(.2*.4^3-.18*.38^3)
Qa = .195*.01*.190
Qb = 2*.1*.01*.2+.195*.18*.01

V = 300*1000

qa = .5*(V*Qa/I)
qb = .5*(V*Qb/I)

%% 2.4

yTild = [5,5,35,35,55];
yTild = yTild/1000;

A = 10*[40,40,50,50,40];
A = A/(1000^2);

ybar = sum(yTild.*A)/sum(A)

B = [40,40,10,10,40]/1000
h = [10,10,50,50,10]/1000

I = (1/12)*B.*h.^3
Ibar = sum(I)+sum(A.*(ybar-yTild).^2)

%% 2.6
clear all

b = 5;%mm
A = 200;
C = 100;
B = 300;

b = b/1000;%m
A = A/1000;
C = C/1000;
B = B/1000;



h = [2*A,2*C];

I = (1/12)*b.*h.^3;

syms x

x = vpa(solve(I(1)*x == I(2)*(B-x)))*1000



