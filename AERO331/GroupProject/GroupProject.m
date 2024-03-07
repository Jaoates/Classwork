clear
close all
clc
syms V r t P L 

assume(V>0)
assume(r>0)
assume(t>0)
assume(P>0)
assume(L>0)
assumeAlso(t<r)


N = 40098.375
L = 1.1
r = .4
V =N*cosd(30)
P =N*sind(30)


r1 = r-t
A = pi*(r^2-r1^2)
I = (pi/4)*(r^4-r1^4)

%% middle

sigxy = (4*V/(3*A))*(r^2+r1*r+r1^2)/(r1^2+r^2)
sigxx = (P/A) 

sig = [[P/A sigxy 0];[sigxy 0 0];[0 0 0]]
% sig = [[r t 0];[0 0 0];[0 0 0]]

[vec,D] = eig(sig)
sige = sqrt(.5*((0-D(2,2))^2+(D(2,2)-D(3,3))^2+(D(3,3)-0)^2))

sigy = 800e6

% solve(sige==sigy)
figure
sige = matlabFunction(sige)
fplot(sige,[.000001,.01])
gca().set("XScale","log")
gca().set("YScale","log")
yline(sigy)
title("middle von mises")

figure
sigxy = matlabFunction(sigxy)
fplot(sigxy,[.000001,.01])
gca().set("XScale","log")
gca().set("YScale","log")
yline(sigy)
title("middle tresca")



%% Top
sigy = 800e6
sigxx = (P/A) + V*L/I
sol = solve(sigxx == sigy)
vpa(sol)

%% 
tin = 0.0019478764268436325198868161624426*1.5;
E = 205e9;
DB = V*L^3/(3*E*I);
DB = subs(DB,t,tin);
DB = vpa(DB)




