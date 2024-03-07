% joshua oates final part 2
clear all
close all
clc

sympref('FloatingPointOutput',true)

CosMat =@(theta) [
    cos(theta)^2            cos(theta)*sin(theta)   -cos(theta)^2           -cos(theta)*sin(theta)
    cos(theta)*sin(theta)   sin(theta)^2            -cos(theta)*sin(theta)  -sin(theta)^2
    -cos(theta)^2           -cos(theta)*sin(theta)   cos(theta)^2           cos(theta)*sin(theta)
    -cos(theta)*sin(theta)  -sin(theta)^2            cos(theta)*sin(theta)  sin(theta)^2
    ];


clear all
close all
clc

sympref('FloatingPointOutput',true)

CosMat =@(theta) [
    cos(theta)^2            cos(theta)*sin(theta)   -cos(theta)^2           -cos(theta)*sin(theta)
    cos(theta)*sin(theta)   sin(theta)^2            -cos(theta)*sin(theta)  -sin(theta)^2
    -cos(theta)^2           -cos(theta)*sin(theta)   cos(theta)^2           cos(theta)*sin(theta)
    -cos(theta)*sin(theta)  -sin(theta)^2            cos(theta)*sin(theta)  sin(theta)^2
    ];

l = @(theta) [
    cos(theta) sin(theta)
    -sin(theta) cos(theta)
    ];

ll =@(theta) [
    l(theta),zeros(2)
    zeros(2),l(theta)
    ];

E = 30e6;
d = .25;
A = pi*(d/2)^2;

EA_L1 = E*A/norm([8 12]);
theta1 = atan(8/12);
K1 = (EA_L1)*CosMat(theta1);

EA_L2 = E*A/8;
theta2 = -pi/2;
K2 = (EA_L2)*CosMat(theta2);

K1 = [
    K1,zeros(4,2)
    zeros(2,6)
    ];

K2 = [
    zeros(2,6)
    zeros(4,2),K2
    ];

K = K1+K2;

syms X F [6 1]

F(3) = 50;
F(4) = 0;
X(1:2) = 0;
X(5:6) = 0;

eqn = F == K*X;
sol = solve(eqn);
F = double(subs(F,sol))
X = double(subs(X,sol));
Xmm = X*10

% Xloc = ll(theta2)*X(1:4)






