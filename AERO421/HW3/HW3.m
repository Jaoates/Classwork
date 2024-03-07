
clear all
close all
clc

addpath("C:\joshFunctionsMatlab\")

mu = 398600;

% syms w(t) r [3,1]
% dw = diff(w);


% current unkowns
syms l mpt 

% cuboid dims
m = 100;%kg
d = 1;%m

% Moment of inertia matrix for the cube
% MOI of cuboid is a lookup formula
I1 = eye(3)*(1/12)*m*(d^2+d^2);

% MOI of one "lump" mass
Ilump = [
    (sind(60)*l)^2 0 0
    0 l^2 0
    0 0 (cosd(60)*l)^2
    ]*mpt;

% total s/c MOI
I = I1+4.*Ilump


% setting up ODE that governs s/c movment
syms theta(t) 

assume(l>0)
assumeAlso(t>0)

% w0 is omega of LVLH to ECI
w0 = 2*pi/(6*3600);

% Tdy as given in the problem statement
T = 1e-5;

%plug in m=2kg
I = subs(I,mpt,2);
Ix = I(1,1);
Iy = I(2,2);
Iz = I(3,3);

%prep for ODE
alpha = 3*w0^2*(Ix-Iz)/Iy;
dtheta = diff(theta);
ddtheta = diff(dtheta);

% from text: ddtheta+alpha*theta = 0
% substitute for: ddtheta+alpha*theta = Tdy/Iz to account for radiation
% pressure torque. Following the derivation steps in the book, it can be 
% seen that putting Tdy/Iz on the left side is sensible.
ODE = ddtheta + alpha*theta - T/Iy == 0;

%initial conditions given in problem statement
icond = [
    theta(0) == 0
    dtheta(0) == 0] ;

%solve ODE in terms of time and l 
theta = dsolve(ODE,icond);
theta = rewrite(theta,"sincos");
disp("theta=")
disp(vpa(simplify(theta),3))

%solve dtheta
dtheta = diff(theta);

%find where dtheta is flat, this is an extreme point, ie we want these
%locations to have an absolute value of less than 5 deg
extrema = solve(dtheta == 0,t);

%plug in the extreme point so that theta is a function of only l
theta = subs(theta,t,extrema);

%solve l such that theta = 5deg as the minimum length of the arm
l = vpa(solve(theta==deg2rad(5)))



