%% Joshua Oates - Exam 2 - Question 5 6 7 11 12
clear all
close all
clc

z = 300;
s = 100;
orientation = 1; % tumbling
year = 2023; % year
d = logspace(-6,3);
flux = zeros(1,length(d));

% inc = 52;
psi = 1.030; % from table in notes at 52deg
for i = 1:length(d)
    flux(i) = Fod(orientation,d(i),z,s,year,psi);
end

figure
hold on
plot(d,flux)
gca().set('xscale','log')
gca().set('yscale','log')
% legend("Background Flux","Actual Flux")
xlabel("Diameter [cm]")
ylabel("Flux [1/(yr*m^2)]")
title("Debris Flux vs Diameter in 2023")

disp("The flux for diameter = 10^-3cm and greater in 2023 is: "+string(Fod(orientation,10^-3,z,s,year,psi))+" [1/(yr*m^2)].")
disp("The flux for diameter = 10^-3cm and greater in 2030 is: "+string(Fod(orientation,10^-3,z,s,2030,psi))+" [1/(yr*m^2)].")
disp("This plot makes sense seeing as how there are significanly less large objects than small ones. I beleive that the bump in the plot from about 10 to 300 cm is from newly dead satalites which would be roughly in this range. They haven't broken up into smaller peices yet. Speaking to the other calculations, it also makes sense that there will be more debris in 2030 than there is in 2023. This is because the model takes into account the increase in new objects from new launches as well as the break up of objects over time in to more numerous smaller objects.")
%% for space shuttle question
A = 735;
t = 1;
F = Fod(1,.8,z,s,year,psi);

PNCF = exp(-F*A*t);
disp("The PNCF for the space shuttle would be: "+string(PNCF)+".")

F = Fod(1,.3,z,s,year,psi);

PNP = exp(-F*A*t);
disp("The PNP for the space shuttle would be: "+string(PNP)+".")




%% functions

function flux = Fod(orientation,d,z,s,t,psi)
% psi from look-up table on average inclination of orbital debris
flux = k(orientation)*H(d)*phi(z,s)*psi*(F1(d)*g1(t)+F2(d)*g2(t));
end

function out = H(d)
% first form from Spenvis, appears to have same output
%     out = 10^(.5*exp(-(log10(d)-0.78)^2/0.406));
out = (10^(exp(-((log10(d)-0.78)/.637)^2)))^.5;
end

function out = phi(z,s)
out = phi1(z,s)/(1+phi1(z,s));
end

function out = phi1(z,s)
out = 10^(z/200 - s/140 - 1.5);
end

function out = F1(d)
out = 1.22e-5*d^-2.5;
end

function out = F2(d)
out = 8.1e10*(d+700)^-6;
end

function out = g1(t)
q1 = 0.02;
q2 = 0.04;

if t<2011
    out = (1+q1)^(t-1988);
else
    out = (1+q1)^23*(1+q2)^(t-2011);
end
end

function out = g2(t)
p = 0.05;
out = 1+p*(t-1988);
end

function out = k(orientation)
arguments
    orientation = 1
end
% case 1 is tumbleling/default
% case 2 is for wake
% case 3 is for ram
switch orientation
    case 1
        out = 1;
    case 2
        out = 0;
    case 3
        out = 3.5;
end
end