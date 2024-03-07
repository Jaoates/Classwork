gam = 1.4
R = 287

X = sqrt(gam/R)

Y = (gam-1)/2

Z = -((gam+1)/(2*gam-2))

De = .42;
Ae = pi*(De/2)^2;
A3 = Ae;
% MFP = mdot*sqrt(Tt3)/(A3*Pt3);
MFP = .0102

syms M
eqn = MFP/X == M + Y*M^7
solve(eqn)
