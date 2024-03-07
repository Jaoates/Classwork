%% 446 final
clear all
close all
clc

%% part 1
re = 6378;%km
z = 700;
r = re+z;
FOR = 2*asind(re/r)

%% part 2
mu = 398600;%km^3/s^2
T = 2*pi/sqrt(mu)*r^(3/2);
T = T/60
Teclipse = T*(FOR/360)

%% part 3

Nd = 300e3/100
FOV = 2*atand(150/700)
lpix = 50;%um
lpix = lpix*1e-6;%m
ldet = lpix*Nd

%% part 4
Ed = 2e-6;% J
v = 6500;%  m/s
pixW = 100;%m
rowRate = v/pixW; %rows/s
Erow = Nd*Ed; %J/row
P0 = Erow*rowRate*4
P20 = P0*1.2
etacool = 1/10;
Pcool = P20/etacool
etasupply = .8;
Psupply = (P20+Pcool)/etasupply




% %% part 5
% % SNR = NEdT
% SNR = (Ad*Df)^.5/(D*A*Ohm*Dv*tau)
% % NEdT =  λ*T2 * (1-exp(-0.0143883*Δλ/λ2))/(SNR*0.0143883)
% 
% NEdT = lam*T^2*(1-exp(-0.0143883*Dlam/lam^2))/(SNR*0.0143883)


%% part 6

sig = 5.6704e-8; % steffan boltzmann
P = (Psupply+150)*1.2
% P = 208.1
epsilon = .8
T1 = 220
T2 = 293
A = P/(epsilon*sig*(T2^4-T1^4))
l = sqrt(A/4)


%% part 7 
P = P/1.2
P = (P +100)*1.2
Tsun = T-Teclipse
Pp = P*(T/Tsun)
eta0 = 245
eta10 = eta0*(1-.02)^10
Ap = Pp/eta10

Ebat = P*Teclipse*60*0.000277778*(1/.5)



%% part 8 
P_CS10 = .6; % proabablility of component success after 10 years
P_MS10 = P_CS10^2; % proababililty of series module success after 10 years. If either component fails, the module will fail
P_MF10 = 1-P_MS10;
n = 4; % number of parallel modules
P_SF10 = P_MF10^n; % probabililty of system failure after 10 years
P_SS10 = 1-P_SF10 % probability of system success after 10 years

disp("You need 4 modules all made from 2 components in series to acheive a proabability of success that is greater than .8")

%% part 9
bit_d = 14;
d0 = Nd*rowRate*2*bit_d
d2 = d0*1.02
dc = d2/1.8
do = dc*T*60

GBo = do*1.25e-10


%% part 10
Xband = 8e9; % Hz
BW = 100e6; % Hz
DR_dl = 400e6; % bps down link

t = do/DR_dl

gt = v*t; %m
gt = gt/1000
FOV = 2*atand((gt/2)/z)
r = z/cosd(FOV/2)

c = 299792458;
lam = c/Xband

BW = c/BW

d = 70*lam/BW

A = pi*(d/2)^2
eta = .8
G = 4*pi*eta*A/lam^2
G = 10*log10(G)

P = 10
EIRP = 10*log10(P)+G

Losses = .3+1+.5+2+.0093*r
CNR = EIRP+31 - Losses
Lfs = 10*log10((lam/(4*pi*r))^2)
A = G*lam^2/(pi*eta*4)
d = 2*sqrt(A/pi)



