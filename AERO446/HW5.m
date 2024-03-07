%% HW 5 - Joshua Oates
clear all
close all
clc

%% Problem 1
clear
bits = 8;
uintMax = 2^bits-1;
myString = '';
for i = 1:bits
    myString = [myString , '1'];
end

% uintMaxCheck = bin2dec(myString)

disp("assuming that we are encoding an unsigned integer (uint), 8 bits can encode the integers 0 to "+string(uintMax))

% Problem 2
Vmax = 1;
lsb = Vmax/(uintMax+1);
disp("The lsb will have a voltage of "+string(lsb)+" V")

% Problem 3
% DR = 20*log(Vmax/lsb)
DR = 20*log(2^bits)

% Problem 4
syms V1
eqn = -6 == 20*log10(V1/Vmax);
V1 = double(solve(eqn,V1))

% Problem 5
T = 290;
kb = 1.38064852e-23;
Gdb = 40;%dB
NF = 6;%dB
F = 10^(NF/10);
G = 10^(Gdb/10);
B = 100e6;

Pnoise = kb*T*B*G*F

% Problem 6
R = 1e3;
syms Inoise Vnoise
eqn = Pnoise == Inoise^2*R;
eqn = [Vnoise == Inoise*R;eqn];
sol = solve(eqn);
Inoise = max(double(sol.Inoise));
Vnoise = max(double(sol.Vnoise))
disp("Since Vnoise > Vlsb, the noise could toggle the lsb")

% Problem 7
P = lsb^2/R

% Problem 8

bits = 9;
uintMax = 2^bits-1;
lsb = Vmax/(uintMax+1)
P2 = lsb^2/R

% Problem 9
Pdiff = 10*log10(P2/P)

% Problem 10
disp("The 9 bit ADC allows for smaller lsb voltage which is an advantage because the ADC will be able to distinguish smaller values from eachother. It also has a thereshold low enough that the noise will not be registered on it.")




%% Functions
% function dB = dec2dB (dec)
% dB = 20*log10(dec);
% end
% function dec = dB2dec(dB)
% dec = 20.^(dB./10);
% end



