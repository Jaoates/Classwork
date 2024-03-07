clear all;close all;clc
addpath("C:\joshFunctionsMatlab\")

%% Problem 1
syms R h m Tx Ty
syms w(t) [3 1]

assume(R,'real')
assume(h,'real')
assume(m,'real')
assume(Tx,'real')
assume(Ty,'real')
assume(w(t),'real')

Ix = (1/12)*m*(3*R^2+h^2);
Iz = m*R^2/2;
I = [[Ix,0,0];[0,Ix,0];[0,0,Iz]];
wx = joshCross(w);
T = [Tx;Ty;0];
ode = diff(w,t) == -inv(I)*(wx*I*w-T);
% cond = w(0) == [0;0;0];
[w1,w2,w3] = dsolve(ode);%,cond)

w1 = simplify(expand(w1))
w2 = simplify(expand(w2))
w3 = simplify(expand(w3))

% clear w1 w2 w3 ode

% a = (Iz-Ix)*w1/Ix;
% ode1 = diff(w1,t,2) == -a*((Ty/Ix)-a*w1);
% 
% [wsol1] = dsolve(ode1)
% [wsol2] = dsolve(ode2)

%% derivative
syms C1 C2 C3 C4 t

diff(C1*cos(C2*t) + C3*cos(C2*t) +C4,t)



%% Problem 2

syms la
syms J [3 3]
syms x cx [3 1]

% cx = conj(x);
J(2,1) = J(1,2);
J(3,1) = J(1,3);
J(3,2) = J(2,3);

% RHS = x'*J*x
% RHS = subs(RHS,conj(x),cx)
% % simplify(x'*J*x)
% 
% x'*x;
% 
% x'*la*x;
% 
% P1 = J(:,1);
% P2 = J(:,2);
% P3 = J(:,3);
% 
% isAlways([cx.'*P1,cx.'*P2,cx.'*P3] == cx.'*J)
% isAlways([cx.'*P1+cx.'*P2+cx.'*P3] == cx.'*J*x)


assume(J,'real')
assume(x'*J*x>0)

[V,Lam]=eig(J)
v1 = simplify(V(:,1))
v2 = simplify(V(:,2))
v3 = simplify(V(:,3))

Lam = simplify(Lam)

% simplify(v1.'*v2)
% simplify(v1.'*v3)
% simplify(v2.'*v3)

%% verify
clear all
J = [...
    [2 -1 0];...
    [-1 3 0];...
    [0 0 1]];

[E,Lam] = eig(J);

disp


