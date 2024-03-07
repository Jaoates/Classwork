% quiz1 - Joshua Oates 
%% Problem 1
% What is the charge per length deposited in Cadmium Telluride if the linear energy transfer per mass
clear

LET = 4; %MeV cm^2 mg^-1
LET = LET*1e6; %eV cm^2 mg^-1
LET = LET*1e3; %eV cm^2 g^-1
rho = 6.2;%g/cm^3
E = 4.43;%eV
q = 1.60217663e-19;% coulombs

Q_x = LET*rho*q/E
answer = Q_x*1e10 % canvas order of magnitude change



%% Problem 2
% What is the absorbance at 2.3 mm if the attenuation length of the material is 10.0 mm?
clear
x = .9;
X = 6.7;

tau = exp(-x/X);
A = -log10(tau)

%% Problem 3
% An astronaut in LEO is exposed to 1.5 Gy of proton radiation and 2.5 Gy of electron radiation equally through the skin, stomach, lungs, and brain. What is the effective dose in Sv to two decimal places?
clear
%[proton;electron]
DT = [2;2.8];
WR = [2;1];
%[skin;stomach;lungs;brain]
WT = [.01;.12;.12;.01];

HT = sum(DT.*WR);

E = (1/length(WT))*sum(WT*HT)


%% Problem 4
% What percentage of 10MeV photons will pass through 2.2 cm silver? Use a linear attenuation of 0.8 cm^-1.
clear

E = 10;% MeV
E = E*1e6;% eV
x = 2.2;%cm
mu = .8;% 1/cm
Ix = E * exp(-mu*x);
percent = Ix/E * 100


%% Problem 5
% What is the exitance of a gray body at 785 K with an emittance of 0.36? 
clear
T = 343;%K
ep = .28;% emittance
sig = 5.6704e-8;
M = ep*sig*T^4



%% Problem 18
% What is the Larmor radius for a hydrogen ion (single proton) moving with a perpendicular velocity of 3.3x10^3 m/s through a magnetic field with average strength 0.5 micro-Tesla?
clear 
m = 9.1093837e-31; %kg (electron)
m = 1.67262192e-27; %kg (proton)
Vt = 5.2*10^3; %m/s
q = 1.60217663e-19;% coulombs
B = 0.8;% micro-Tesla
B = B*1e-6;%tesla
rl = m*Vt/(q*B)
