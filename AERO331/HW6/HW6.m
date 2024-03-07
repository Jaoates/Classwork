%% HW6 - A331 - Joshua Oates
clear all
close all
clc

%% Problem 1.1
clear all
syms d G Mx A
assume(d>0)
assume(G>0)
assume(Mx>0)

% cross section 1 is square
% cross section 2 is circle

% A = .01*d^2;
A_bar1 = d^2;
A_bar2 = pi*(d/2)^2;

s1 = 4*d;
s2 = pi*d;

t1 = A/s1;
t2 = A/s2;

q1 = Mx/(2*A_bar1);
q2 = Mx/(2*A_bar2);

theta1 = (1/(2*A_bar1))*(q1/(G*t1))*s1
TR1 = Mx/theta1

theta2 = (1/(2*A_bar2))*(q2/(G*t2))*s2
TR2 = Mx/theta2

isAlways(TR1==TR2)

%% Problem 1.2
clear all
syms d G Mx t
assume(d>0)
assume(G>0)
assume(Mx>0)
assume(t>0)

% t = t1 = t2

s1 = 4*d;
s2 = pi*d;

A1 = t*s1;
A2 = t*s2;

A_bar1 = d^2;
A_bar2 = pi*(d/2)^2;

q1 = Mx/(2*A_bar1);
q2 = Mx/(2*A_bar2);

theta1 = (1/(2*A_bar1))*(q1/(G*t))*s1
TR1 = Mx/theta1;

theta2 = (1/(2*A_bar2))*(q2/(G*t))*s2
TR2 = Mx/theta2;

isAlways(TR1>TR2)



%% Problem2
clear all
syms q1 q2 q3 Mx

assume(q1>0)
assume(q2>0)
assume(q3>0)
assume(Mx>0)

qw1 = q1-q2;
qw2 = q2-q3;

G = 4e6;
t1 = .05;
t2 = .1;

A1 = pi/2;
A2 = 6;
A3 = 3;

s1 = pi;
s2 = sqrt(2^2+3^2);

eqn1 = Mx/2 == q1*A1 + q2*A2 + q3*A3;
eqn2 = (1/A1)*((q1*s1/t1)+(qw1*2/t1)) == (1/A2)*(-qw1*2/t1 + 2*3*q2/t2 + qw2*2/t1);
eqn3 = (1/A1)*((q1*s1/t1)+(qw1*2/t1)) == (1/A3)*(-qw2*2/t1 + q3*(s2+3)/t1);

%a
sol = solve([eqn1,eqn2,eqn3]);

q1 = sol.q1;
q2 = sol.q2;
q3 = sol.q3;

q1 = vpa(q1)
q2 = vpa(q2)
q3 = vpa(q3)


qw1 = q1-q2;
qw2 = q2-q3;

%b
sig1 = q1/t1
sig2 = q2/t2
sig3 = q3/t1

sigw1 = qw1/t1
sigw2 = qw2/t1

%c
theta = vpa( (1/(2*A1*G)) * ((q1*s1/t1)+(qw1*2/t1)) )

%d
TR = Mx/theta

%e
TRclass = 6e6;

TR_ratio = TR/TRclass





