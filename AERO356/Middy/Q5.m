%% Question 5
clear all
close all
clc
% part 1

q = 1.60217663e-19;% coulombs
m = 9.1093837e-31; %kg (electron)
ke = 100;%keV
ke = ke*1e3; % eV
ke = ke*q;% J

Vt = sqrt(2*ke/m);

B = .15e-4;%T
rl = m*Vt/(q*B)

% part 2

m = 1.67262192e-27; %kg (proton)
Re = 6378.1e3;%m

c = 299792458;

ke = B*c*10*Re %eV