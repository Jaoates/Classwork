% 446 Hw 4 - Joshua Oates
clear all
close all
clc
% Problem 1
%{
Calculate the gain (in dB) of a 1-meter diameter parabolic antenna operating at a
frequency of 300 MHz. Assume an antenna efficiency of 60%
%}
clear

c = 299792458;
f = 300e6;
lam = c/f;
eta = .6;
D = 1;

Aeff = pi*(D/2)^2;
G = eta*4*pi*Aeff/lam^2;
GdB = dec2dB(G);
disp("Problem 1")
disp("f = 300MHz:") 
disp("at the nominal wavelength, the gain in dB is: "+string(GdB))

% Problem 2
%{
Consider the above antenna operated at 10% bandwidth, that is from -5% to +5% of
the nominal frequency. Is the gain symmetric around the nominal 300 MHz? Why or
why not?
%}
f = [.95*f,f,1.05*f];
lam = c./f;
G = eta*4*pi*Aeff./lam.^2;
GdB = dec2dB(G);
dGdB = abs(diff(GdB));
disp("Problem 2")
disp("The difference from the nominal gain the +- 5% gain in the positive direction is: "+string(dGdB(2)))
disp("The difference from the nominal gain the +- 5% gain in the negative direction is: "+string(dGdB(1)))
disp("This means that the antenna is more effective at +5% wavelength than it is at -5%")
disp("The behavior is asymetric around the nominal value")

% Problem 3
%{
Repeat problem 1 at 1 GHz and 10 GHz.
%}
% Gdb = 10*log10(Pi/Po)
c = 299792458;
f = 1e9;
lam = c/f;
eta = .6;
D = 1;

Aeff = pi*(D/2)^2;
G = eta*4*pi*Aeff/lam^2;
GdB = dec2dB(G);
disp("Problem 3")
disp("f = 1GHz:")
disp("at the nominal wavelength, the gain in dB is: "+string(GdB))


f = 10e9;
lam = c/f;
eta = .6;
D = 1;

Aeff = pi*(D/2)^2;
G = eta*4*pi*Aeff/lam^2;
GdB = dec2dB(G);
disp("f = 10GHz:")
disp("at the nominal wavelength, the gain in dB is: "+string(GdB))


% Problem 4
%{
Convert 0.1 dB to a percentage
%}
leveldB = .1;
level = dB2dec(leveldB);
level = level*100;
disp("Problem 4")
disp(".1dB in percent is "+ string(level)+" %")

% %% Problem 5
% Problem 6
clear 
GA = 30;%dB
GB = 20;%dB
GC = 13;%dB
G0dB = [GA GB GC];
GtotdB = sum(G0dB);
G0 = dB2dec(G0dB);
% FTOTAL= F1 + [F2 - 1]/G1 + [F3 - 1]/[G1G2]
F0 = [2,1.6,1.4];
ord = [
1 2 3
1 3 2
2 1 3
2 3 1
3 1 2
3 2 1
];

for i = 1:length(ord)
    Gm(i,:) = [G0(ord(i,1)),G0(ord(i,2)),G0(ord(i,3))];
    Fm(i,:) = [F0(ord(i,1)),F0(ord(i,2)),F0(ord(i,3))];
end


Ftot = @(F,G) F(1) + (F(2)-1)/G(1) + (F(3)-1)/(G(1)*G(2));

for i = 1:length(Gm)
    Ftotm(i) = Ftot(Fm(i,:),Gm(i,:));
end
Ftotm =Ftotm';
[Ftotmin,I] = min(Ftotm);
NFtotmin = dec2dB(Ftotmin);
NFtotm = dec2dB(Ftotm);
disp("Problem 6")
disp("The total receiver gain in dB is: "+string(GtotdB))
disp("The minimimum noise configuration is: "+sprintf('%d %d %d',ord(I,:)))
disp("NF in dB for this configuration is: "+string(NFtotmin))

% Problem 7
i = 3;% 3rd permutation is 213 is BAC
syms FC 
F = [Fm(i,1),Fm(i,2),FC];
G = Gm(i,:);

Ftot = Ftot(F,G);
eqn = Ftot == Ftotm(i)*dB2dec(0.1);
FC = double(solve(eqn,FC));
NFC = dec2dB(FC);
disp("Problem 7")
disp("NFC could go as high as "+string(NFC)+" dB and the overall noise increase would be less than .1 dB")


%Problem 8
dT = (F0-1)*290';
disp("Problem 8")
disp("The equivilent delta temperature for all three gains in order ABC in [K] is as follows: ")
disp(dT)

% %% Problem 9
% disp("Problem 9")



% Functions
function dB = dec2dB (dec)
dB = 10*log10(dec);
end
function dec = dB2dec(dB)
dec = 10.^(dB./10);
end