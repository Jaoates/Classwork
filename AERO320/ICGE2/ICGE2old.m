clear all
close all
clc

addpath("C:\joshFunctionsMatlab\")


%% 1
L = 4; %m
R = .25; %m
h = 1; %m
mb = 600; %kg
mn = 50; %kg
mw = 10; %kg
r = .2; %m
t = .04; %m

Vb = pi*R^2*L;
sigb = mb/Vb;

Vn = pi*R^2*h/3;
sign = mn/Vn;

c = ((.5*L)*mb+(L+.25*h)*mn)/(mn+mb);

% adapted from the previous t:

% J of cylinder from bottom center face
Jcyl =@ (R0,h0,m0) [...
[(m0*(3*R0^2 + 4*h0^2))/12,                      0,         0];...
[                     0, (m0*(3*R0^2 + 4*h0^2))/12,         0];...
[                     0,                      0, (R0^2*m0)/2]];

    
% J of cone from tip
Jcon =@(R0,h0,m0)[ ...
[(3*m0*(R0^2 + 4*h0^2))/20,                             0,              0];...
[                            0, (3*m0*(R0^2 + 4*h0^2))/20,              0];...
[                            0,                       0, (m0*3*R0^2)/10]];

I =@(J,m,rc) J + m*joshCross(rc)*joshCross(rc);

Jw = Jcyl(r,t,mw);
Jb = Jcyl(R,L,mb);
Jn = Jcon(R,h,mn);

Iw = I(Jw,mw,[0;0;.5*t]);
Ib = I(Jb,mb,[0;0;c]);
In = I(Jn,mn,[0;0;-(L-c)-h]);
Ir = Ib+In

%% Hammer time


r21cyl = [0 0 -L/2];
r21cylcross = joshCross(r21cyl);
Ib = Jb + (mb*r21cylcross*r21cylcross);
% J matrix from cone tip

Jcone = @(r,h,m) [(3*m*((r^2)+(4*(h^2))))/20 0 0;...
    0 (3*m*((r^2)+(4*(h^2))))/20 0;...
    0 0 (m*3*(r^2))/10];
Jn = Jcone(R,h,mn);
r21cone = [0 0 0.75*(h)];
r21tipcross = joshCross(r21cone);
In = Jn + (mn*r21tipcross*r21tipcross);

% I matrix about rocket
r21nose = [0; 0; -((L+.25*h)-c)];
r21nosecross = joshCross(r21nose);

Ir_nose = In - (mn*r21nosecross*r21nosecross);

r21body = [0 0 c-2];
r21bodycross = joshCross(r21body);

Ir_body = Ib - (mb*r21bodycross*r21bodycross);
Ir2 = Ir_nose + Ir_body

% I matrix about reaction wheel
Icyl = @(m,r,h) [(m*((3*(r^2))+(h^2)))/12 0 0;...
    0 (m*((3*(r^2))+(h^2)))/12 0;...
    0 0 (m*(r^2))/2];
Iw2 = Icyl(mw,r,t);




% % Center of Mass
% L = 4;
% R = 0.25;
% h = 1;
% mw = 10;
% t = 0.04;
% r = 0.2;
% mn = 50;
% mb = 600;
% mtot = mn + mb;
% (L+.25*h) = 4.25;
% 2 = 2;
% ztot = (L+.25*h) + 2;
% Zm_cone = 1200;
% Zm_body = 212.5;
% Zm_tot = Zm_cone + Zm_body;
% c = Zm_tot/mtot;    % From bottom
% % Adapted matrices from recent t
% % J matrix from cylinder bottom

% Jcyl = @(r,h,m) [(m*(3*(r^2)+(4*(h^2))))/12 0 0;...
%     0 (m*(3*(r^2)+(4*(h^2))))/12 0;...
%     0 0 (m*(r^2))/2];

% J = Jcyl(R,L,mb);
%%
clear all
syms hc w wrel [3 1]
syms Ir Iw [3 3]

wdot = diff(w);
wreldot = diff(wrel);
hcdot = diff(hc);

wx = joshCross(w);
wrel = joshCross(wrel);

Td = Ir*wdot+wx*Ir*w + Iw*wreldot+wrel*Iw*wrel;
Td = simplify(Td);

Td = subs(Td,w1,0);
Td = subs(Td,w2,0);
Td = subs(Td,wrel1,0);
Td = subs(Td,wrel2,0);

Td = simplify(Td);
LHS = Td;
clear Td
syms Td [3 1]
eq = sym(Td == LHS);
eq = solve(eq,wdot)