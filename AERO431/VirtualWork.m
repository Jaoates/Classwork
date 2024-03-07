%% unit 4 section 2
% Joshua Oates
clear all
close all
clc

%% 14.13
clear
%{
          |
          V
<--  -----
    |    /|
    |  /  |
    |/    |
    -------
    O      *
%}

% virtual load
P2 = 0;%kN
P1 = -1;%kN
a = 2;%m
b = 1.5;%m

n1 = 0;
n2 = P2;
theta = atan(b/a);
n3 = -P2/cos(theta);
n4 = P2;
n5 = -P1-sin(theta)*n3;

% real load
P2 = 4e3;%kN
P1 = 5e3;%kN
a = 2;%m
b = 1.5;%m

N1 = 0;
N2 = P2;
theta = atan(b/a);
N3 = -P2/cos(theta);
N4 = P2;
N5 = -P1-sin(theta)*N3;

A = 400;%mm
A = A*(1e-3)^2;
E = 200e9;
delta = n5*N5*b/(A*E);

%% 14.15
clear

% virtual load
P2 = 1;%kN
P1 = 0;%kN
a = 2;%m
b = 1.5;%m

n1 = 0;
n2 = P2;
theta = atan(b/a);
n3 = -P2/cos(theta);
n4 = P2;
n5 = -P1-sin(theta)*n3;

% real load
P2 = 10e3;%kN
P1 = -5e3;%kN
a = 2;%m
b = 1.5;%m

N1 = 0;
N2 = P2;
theta = atan(b/a);
N3 = -P2/cos(theta);
N4 = P2;
N5 = -P1-sin(theta)*N3;

L1 = b;
L2 = a;
L3 = sqrt(a^2+b^2);
L4 = a;
L5 = b;

A = 400;%mm
A = A*(1e-3)^2;
E = 200e9;
delta = n1*N1*L1+n2*N2*L2+n3*N3*L3+n4*N4*L4+n5*N5*L5;
delta = delta/(E*A)




%% 14.16
clear
% virtual load
P1 = 1;
P2 = 0;
a = 1.5;
b = 2;

n1 = P1+P2;
n2 = (P1*2*a+P2*a)/b;
theta = atan(b/a);
n3 = -(P1+P2)/sin(theta);
% n5 = -P2+(sin(theta)*n3)/sin(theta);
n5 = P1/sin(theta);
n4 = -cos(theta)*n5;

% real load
P1 = 20;
P2 = 20;
a = 1.5;
b = 2;

N1 = P1+P2;
N2 = (P1*2*a+P2*a)/b;
theta = atan(b/a);
N3 = -(P1+P2)/sin(theta);
% N5 = -P2+(sin(theta)*N3)/sin(theta);
N5 = P1/sin(theta);
N4 = -cos(theta)*N5;

L1 = b;
L2 = a;
L3 = sqrt(a^2+b^2);
L4 = 2*a;
L5 = L3;

d = 30;%mm
A = pi*(d/2)^2;%mm
A = A*(1e-3)^2;
E = 200e9;

delta = n1*N1*L1+n2*N2*L2+n3*N3*L3+n4*N4*L4+n5*N5*L5;
delta = delta/(E*A)

%% 14.17
E = 13e9;
A = 180e-3*120e-3




