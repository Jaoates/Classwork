% joshua oates final part 1
clear all
close all
clc

sympref('FloatingPointOutput',true)


L = 600e-3;
E = 50e9;
EA_L = [2400 1800 1200]*1e-6*E/L;
P = 28e3;

k = [1 -1; -1 1];


K1 = EA_L(1)*k
K2 = EA_L(2)*k
K3 = EA_L(3)*k

K = zeros(4);
K(1:2,1:2) = K(1:2,1:2)+ K1;
K(2:3,2:3) = K(2:3,2:3)+ K2;
K(3:4,3:4) = K(3:4,3:4)+ K3



syms F X [4 1]
F(4) = P
F(2:3) = 0
X(1) = 0

eqn = F == K*X;
sol = solve(eqn);
F = double(subs(F,sol))
X = double(subs(X,sol));
Xmm = X*1000

B = [-1 1]*(1/L);
epsilon(1) = B*X(1:2);
epsilon(2) = B*X(2:3);
epsilon(3) = B*X(3:4)

sig = E*epsilon








