clear all
close all
clc

addpath("C:\joshFunctionsMatlab\")

% this file will calculate the initial orientation of the sc

% inertial basis vectors expressed in inertial frame (ECI)
xi = [1;0;0];
yi = [0;1;0];
zi = [0;0;1];

% Rotation matrix given with assignment
Ci0 =[
    0   -0.1466    0.9892
    0    0.9892    0.1466
    -1.0000         0         0];

% LVLH basis vectors
x0 = Ci0*xi;
y0 = Ci0*yi;
z0 = Ci0*zi;

% get the quaternion
[eta,eps] = joshRotM2Quat(Ci0);

% get euler angles
[t3,t2,t1] = quat2angle([eta;eps]');


mu = 398600.4418;
h=53335.2;
ecc=0;
raan = 0;
inc = deg2rad(98.43);
aop = 0;
theta = 0;
a = (h^2/mu)/(1-ecc);

[r0,v0]  = joshCOE2rv(a,ecc,theta,inc,raan,aop,mu);



% output all values
initial.t1i0 = t1;
initial.t2i0 = t2;
initial.t3i0 = t3;


initial.t1i0 = 1.5708;
initial.t2i0 = -1.4237;
initial.t3i0 = -1.5708;

initial.etai0 = eta;
initial.epsi0 = eps;
initial.Ci0 = Ci0;
initial.v0 = v0';
initial.r0 = r0';
initial

% save values
save("initialOrientation.mat","initial","-mat")



