clear all
close all
clc

addpath("C:\joshFunctionsMatlab\")

% Problems 
% 8.2
% 8.4
% 8.6
% 8.7
% 8.12
% % 8.16  


%% 8.2
mu_s = 132712440018;
r1 = 227.9e6;
r2 = 778.5e6;
v1 = sqrt(mu_s/r1);
v2 = sqrt(mu_s/r2);

[dv1,dv2,dv,T,ht,ecct,vt1,vt2] = joshHomann(r1,v1,r2,v2,mu_s);

disp("------------P8.2------------")
disp("My calculations have the following results:")
disp("dv: "+string(dv)+" km/s")
disp("H/C: This dv seems reasonable for an interplanetry launch to mars.")
disp("H/C: ecc of transfer: "+string(ecct)+"seems to be reasonable for this trajectory (is ellpitical but not too elliptical).")




%% 8.4
clear all
mu_s = 132712440018;
r1 = 227.9e6;
r2 = 778.5e6;
n1 = sqrt(mu_s/r1^3);
n2 = sqrt(mu_s/r2^3);
Tsyn = 2*pi/(n1-n2);
Tsyn = Tsyn/(60*60*24); %solar days


disp("------------P8.4------------")
disp("My calculations have the following results:")
disp("T-syn Mars/Jupiter: "+string(Tsyn)+" solar days")
disp("Note: solar days was the prefered unit in the book.")
disp("H/C: since this is larger than the orbital period of Mars, and Jupiter is moving much slower than Mars, this T-syn seems reasonable.")
%% 8.6
clear all
mu_sun = 132712440018;
% rj = 778.5e6;
rs = 1.433e9;
ru = 2.872e9;
rn = 4.495e9;
% mu_j = 126686534;
mu_s = 37931187;
mu_u = 5793939;
mu_n = 6836529;

% rsoij = rj*(mu_j/mu_sun)^(2/5);
rsois = rs*(mu_s/mu_sun)^(2/5);
rsoiu = ru*(mu_u/mu_sun)^(2/5);
rsoin = rn*(mu_n/mu_sun)^(2/5);

disp("------------P8.6------------")
disp("My calculations have the following results:")
disp("SOI radii (km):")
disp("Saturn: "+ rsois)
disp("Uranus: "+ rsoiu)
disp("Neptune: "+ rsoin)
disp("H/C: these large SOI values correspond to gas giants.")
disp("H/C: these values are close to published valeus.")

%% 8.7
clear all
mu_e = 398600;
r_e = 6378; % km
rpark = 200+r_e;

mu_sun = 132712440018;
ra = 147.4e6;
rp = 120e6;
v1 = sqrt(mu_sun/ra);
a2 = (ra+rp)/2;
ecc2 = (ra-rp)/(ra+rp);
h2 = sqrt(mu_sun*ra*(1+ecc2*cos(pi)));
va = h2/ra;
vp = h2/rp;
vinf = v1-va;

vpark = sqrt(mu_e/rpark);
vbo = sqrt(vinf^2+2*mu_e/rpark);

dv = vbo-vpark;


disp("------------P8.7------------")
disp("My calculations have the following results:")
disp("Delta-V required: "+string(dv)+" km/s")
disp("Excess-V: "+string(vinf)+" km/s")
disp("H/C: ecc: "+string(ecc2)+" is within a sensible range for this elliptical orbit.")
disp("H/C: Delta-V seems to be reasonable for this manouver.")




%% 8.12
clear all
mu_s = 132712440018;
mu_j = 126686534;
r_j = 778.5e6;
r_e = 149.6e6;
v_j = sqrt(mu_s/r_j);
v_e = sqrt(mu_s/r_e);

rsoij = 4.8215e7;

[dv1,dv2,dv,T,ht,ecct,vt1,vt2] = joshHomann(r_e,v_e,r_j,v_j,mu_s);

vsc1 = [vt2,0];
z = 200000;
rpj = 71490+z;
v_j = [1,0]*v_j;
vinf1 = vsc1-v_j;
sinf = norm(vinf1);
ecc = 1+(rpj*sinf^2)/mu_j;
h = mu_j*sqrt(ecc^2-1)/sinf;
B = acos(1/ecc);
d = 2*B;
vinf2 = sinf*[cos(d),-sin(d)];

vsc2 = vinf2+v_j;
dv = vsc2-vsc1;
dv = norm(dv);

vsc2 = [vsc2,0];
rsc2 = [0,r_j,0];
[a2,ecc2,theta2,inc2,raan2,aop2,h2,T2,E2] = joshCOE(rsc2,vsc2,mu_s,"magnitude");

disp("------------P8.12------------")
disp("My calculations have the following results:")
disp("Delta-V: "+string(dv)+" km/s")
disp("a: "+string(a2)+" km")
disp("ecc: "+string(ecc2))
disp("H/C: the imparted Delta-V is well within the range for a Jupiter flyby, at a high altitude this relatively small (for jupiter) Delta-V makes sense.")
disp("H/C: the magnitude of Excess-V for arrival and departure is the same.")

%% 8.16  
clear all

mu_sun = 132712440018;

tString = "August 15, 2005";
t = datetime(tString);
[~,~,j2000_1]= joshJulian(t);

tString = "March 15, 2006";
t = datetime(tString);
[~,~,j2000_2]= joshJulian(t);

T0_1 = j2000_1/36525; % confirmed T0 is curtis
T0_2 = j2000_2/36525;

% T0_1 = (jd1-2451545)/36525;
% T0_2 = (jd2-2451545)/36525;


% clear t tString j0_2 j0_1

% coe1 = abercrombyAERO351planetary_elements2(3,T0_1);
% 
% a1 = coe1(1);
% ecc1 = coe1(2);
% inc1 = deg2rad(coe1(3));
% raan1 = deg2rad(coe1(4));
% w_hat1 = deg2rad(coe1(5));
% L1 = deg2rad(coe1(6));
% aop1 = deg2rad(w_hat1-raan1);
% Me1 = deg2rad(L1 - w_hat1);
% [theta1,~,E1] = joshAnomalyCalculator(ecc1,Me1,"Me")
% 
% coe2 = abercrombyAERO351planetary_elements2(4,T0_2);


[coe1, r1, vearth, jd1] = curtisPlanet_elements_and_sv(3, 2005, 8, 15, 0, 0, 0);
[coe2, r2, vmars, jd2] = curtisPlanet_elements_and_sv(4, 2006, 3, 15, 0, 0, 0);



coefs = 15;
[Cc,Sc]=joshStumpffCoeffs(coefs);
C = @(z) sum(Cc.*joshStumpffZ(z,coefs));
S = @(z) sum(Sc.*joshStumpffZ(z,coefs));



dt = (j2000_2-j2000_1)*24*3600;

z = cross(r1,r2);
z = z(3);

theta1 = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
theta2 = 2*pi - theta1;

[fz,y,A,z,flag,glag,gdotlag] = joshfLambert(norm(r1),norm(r2),dt,theta1,mu_sun,Cc,Sc);

v1 = (1/glag)*(r2-flag*r1);
v2 = (1/glag)*(gdotlag*r2-r1);

vinf1 = abs(norm(v1) - norm(vearth));
vinf2 = abs(norm(v2) - norm(vmars));
clear A C Cc coe1 coe2 coefs dt flag fz gdotlag glag j2000_1 j2000_2 jd1 jd2 mu_sun r1 r2 S Sc t T0_1 T0_2 theta1 theta2 tString y z

r_e = 6378; % km
mu_e = 398600;
r_m = 3396;
mu_m = 42828;

z = 190;
rbo = z+r_e;
mu_e = 398600;
v0 = sqrt(mu_e/rbo);
vbo = sqrt(vinf1^2+(2*mu_e/rbo));
dv1 = abs(vbo-v0);

zp2 = 300;
rp2 = zp2 + r_m;

P = 35;% hr
P = P*60*60; % sec
a = (P*sqrt(mu_m)/(2*pi))^(2/3);
ra2 = 2*a-rp2;
ecc2 = (ra2-rp2)/(ra2+rp2);

h2 = sqrt(mu_m*rp2*(1+ecc2));
vp2 = h2/rp2;
vba = sqrt(vinf2^2+(2*mu_m/rp2));
dv2 = norm(vba-vp2);

dv = dv1+dv2;

disp("------------P8.16------------")
disp("My calculations have the following results:")
disp("Delta-V: "+string(dv)+" km/s")
disp("H/C: 'long way' and 'short way' are almost the same, this makes sense since arrival is almost on the opposite side of the sun from deprature.")
disp("H/C: Delta-V is in a resonable range for a transfer to Mars.")
disp("Note: My answer doesn't match exactly with the book. In debugging I hard-coded in the book's 'r' vectors and noticed that my lambert's solver is working and also my coes->'r'&'v' vectors function is working. I also double checked my departure and arrival times. I think there may be something wrong with Curtis's 'rates' table or else there is just a decent amount of error in propogating planet's paths this way.")

%% dependancies


disp("------------Dependancies------------")
disp("My code uses the following functions: ")
depends = matlab.codetools.requiredFilesAndProducts('C:\AERO351\A351HW4\HW4');
disp(depends')


