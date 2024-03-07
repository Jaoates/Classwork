clear all
close all
clc
addpath('C:\joshFunctionsMatlab\')

num = sqrt(3)/2

C21 = [[.5 num 0];[0 0 1];[num -.5 0]]
joshIsOnes(joshIsRotM(C21))

x = [1;0;0] % 2_2
y = [0;1;0]
z = [0;0;1]

x1_2 = C21*x % 1_2
y1_2 = C21*y
z1_2 = C21*z

a = dot(x,(2*x1_2-num*z1_2))
b  = dot(y1_2,(x-sqrt(3)*z))

%% 2
clear all

[Cx,Cy,Cz] = joshAxisRotation();
syms a B wt t0 t3

theta2 = pi/2-B
% theta3 = t0 + wt
theta3 = t3

C21 = simplify(Cz(theta3))
C32 = simplify(Cy(theta2)*Cz(a))
C31 = C32*C21

syms Rs
syms rho [3 1]
R3 = [0;0;1]*Rs
r3 = R3+rho

r1 = simplify(C31.'*r3)





