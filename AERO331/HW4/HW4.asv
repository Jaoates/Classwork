%% housekeeping
clear all
close all
clc

%% Problem 2

% sigma vector
syms sig eps [6 1]
syms u [3 1]
syms C11 C12 C44 k1 k2 nu E alpha T 
syms X [3 1]

assumeAlso(X,'real')
assumeAlso(E>0)
assumeAlso(nu>0)
assumeAlso(alpha>0)
assumeAlso(T>0)



% C44 = .5*(C11-C12);

% k1 =(E/((1+nu)*(1-2*nu))); 
C = k1*[
[C11 C12 C12 0   0   0]
[C12 C11 C12 0   0   0]
[C12 C12 C11 0   0   0]
[0 0 0       C44 0   0]
[0 0 0 0         C44 0]
[0 0 0 0 0           C44]
];
C = subs(C,C11,1-nu);
C = subs(C,C12,nu);
C = subs(C,C44,.5*(1-2*nu));

eps(1) = 0;
eps(2) = 0;
eps(6) = 0;

eps(4) = 2*eps(4);
eps(5) = 2*eps(5);

sig = C*eps - k2*[1;1;1;0;0;0];
sig = simplify(sig)

% this is the trace
sig_h = (1/3)*(sig(1)+sig(2)+sig(3))
sig_dev = sig-sig_h*[1;1;1;0;0;0]

sig_e = ((3/2)*sum(sig_dev.^2))^0.5
sig_e = simplify(sig_e)

% symetric
for i = 1:6
    for k = 1:6
        C(k,i) = C(i,k);
    end
end
C = simplify(C);

clear C11 C12 C44

% eqn1 = sig == C*eps - ((E*alpha*T)/(1-2*nu))*[1;1;1;0;0;0];


%%
E = 10^7;
sigy = 50*1000;
a = 10^-5;
T = 100;
nu = 1/3;
k2 = E*a*T/(1-2*nu)
syms p0

nu*(k2/(1-nu))-(p0/(1-nu))-k2


