% practice exam 2
syms k1 nu Eyy Ezz Exy Exz Eyz C11 C12 C44 k2

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

E = [0;Eyy;Ezz;0;0;0];


sig = simplify(expand(C*E))-k2*[1;1;1;0;0;0]
simplify(sig)


%% Problem 1
clear
syms theta k

assume(k>0)
assume(theta>=0)

sig = k*[[0 0 0];[0 0 0];[0 0 1]]
n = [sin(theta);0;cos(theta)]

T = sig*n
Tn = dot(T,n)*n

Tt = T-Tn








