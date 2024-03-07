clear all
close all
clc

addpath("C:\joshFunctionsMatlab\")
%%
syms x m L t;
r = [x;0;0];
rx = joshCross(r);
rxx = rx*rx
I = -int(rxx*m/L,x,[-L/2,L/2])
[Cx,Cy,Cz] = joshAxisRotation();
% J1 = Cz(0+t)*J
% J2 = Cz(pi*(2/3)+t)*J
% J3 = Cz(pi*(4/3)+t)*J
% I = J1+J2+J3
I2 = Cz(t)*I

%%
function dX = odeFun2(t,X)
w = X(1:3);
% wrel = X(4:6);
wrel = zeros(3,1);
E = X(7:9);

%%%%%%%
dwrel =[0;0;0];
%%%%%%%
Td = [0;1;0]*a;
%%%%%%%



mat = [[cos(E(2)),      sin(E(2))*sin(E(1)),            sin(E(2))*cos(E(1))];...
    [0,                 cos(E(2))*cos(E(1)),           -cos(E(2))*sin(E(1))];...
    [0,                 sin(E(1)),                      cos(E(1))           ]]*(1/cos(E(2)));

wx = joshCross(w);
ww=(w+wrel);
wwx = joshCross(ww);
dw = -inv(Ir+Iw)*(wx*Ir*w+Iw*dwrel+wwx*Iw*ww+Td);

Edot = mat*w;

dX = [dw;dwrel;Edot];
end