%% clean up
clear all
close all
clc
%% BLT1



%% SP1
clear all
close all
syms t x y
% U = @(t) [10*cos(10*t)*exp(-t);0];
Ux = 10*cos(10*t)*exp(-t);
Uy = 0;
psiy = int(Ux,y);
psix = -int(0,x);

yfun = matlabFunction(psiy);
yfun2 =@(t) .1.*cos(t.*1.0e+1).*exp(-t).*1.0e+1;
% xfun = matlabFunction(psix);
xfun = @(t) 0;
figure
fplot(xfun,yfun2,[0,10])
title("steamlines X Y")
xlabel("x [m]")
ylabel("y [m]")

figure
fplot(xfun,yfun2,[0,10])
title("pathlines X Y")
xlabel("x [m]")
ylabel("y [m]")

t = @(t) t;
figure
fplot3(xfun,yfun2,t,[0,10])
title("pathlines X Y t")
xlabel("x [m]")
ylabel("y [m]")
zlabel("t [s]")

disp("the flows are unsteady beacuse e pathlines change with respect to time. It looks like it aproaching a steady flow around t = 8s.")

%% EA1
syms t r R Uinf
psi = (r-R^2/r)*Uinf*sin(t);
ur = diff((1/r)*psi,t);
ut = -diff(psi,r);

phir = int(ur,r);
phir = simplify(rewrite(phir,'sincos'));
phit = int(ut*r,t);
isAlways(phir==phit);

phi = phir;
clear phit phir

int(cos(t),t,[0,2*pi])
int(sin(t)^2,t,[0,2*pi])
int(sin(t)^2*cos(t),t)

syms Pinf rho
P = Pinf + .5*rho*Uinf^2-.5*rho*(-2*Uinf*sin(t))^2;
P2 = Pinf+.5*rho*Uinf^2*(1-4*sin(t)^2);
isAlways(P==P2)
eq = -R*cos(t)*P;
int(eq,t,[0,2*pi]);
