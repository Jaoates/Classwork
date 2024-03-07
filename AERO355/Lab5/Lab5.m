clear all
close all
clc
dat = importdata("2017-5-24_AERO353_vibelab.xlsx").Sheet1;
cols = ["X-Data [Hz]","Reference", "Article 1 3749 (M - Filtered) [g]", "Phase [Deg]", "Article 2 8361 (M - Filtered) [g]", "Phase [Deg]","Article 3 8482 (M - Filtered) [g]", "Phase [Deg]","Control [g]"];
dat = dat(:,1:4);
cols = cols(1:4);
u = symunit;
%% a priori easy
Q = 33.3;

E = 29000000*u.psi;
L = 11*u.in;
d = .25*u.in;
r = d/2;
massm = .13*u.lbm;
massr = .09*u.lbm;

I = pi*r^4/4;

keasy =3*E*I/L^3;
Fneasy = (1/(2*pi))*sqrt(keasy/massm);
Fneasy = vpa(simplify(unitConvert(Fneasy,"US")))


%% a priori hard
Q = 33.3;

E = 29000000*u.psi;
L = 11*u.in;
d = .25*u.in;
r = d/2;

weightm = .13*u.lbf;
weightr = .09*u.lbf;
weight = weightr + weightm;


massm = .13*u.lbm;
massr = .09*u.lbm;
mass = massm+massr;

I = pi*r^4/4;

P = weightm;
x = L;
h = L;

w = weightr/L;

delm = ((P*x^2)/(6*E*I))*(3*h-x);
delr = ((w*x^2)/(24*E*I))*(x^2+6*L^2-4*L*x);
del = delm + delr;


khard = weight/del;
Fnhard = (1/(2*pi))*sqrt(khard/mass);
Fnhard = vpa(simplify(unitConvert(Fnhard,"US")))

maxDisplacementEasy = vpa(simplify(unitConvert(Q*del,"US")))



%% Physical testing
Fref = dat(1,2);
[Fmax, iFmax] = max(dat(:,3));
figure 
hold on
plot(dat(:,1),dat(:,3))
plot(dat(iFmax,1),dat(iFmax,3),".")
legend("All","Max")
xlabel("Frequency [Hz]")
ylabel("Amplitude of Acceleration [gee]")

Q = dat(iFmax,3)/Fref
Fn = dat(iFmax,3)
