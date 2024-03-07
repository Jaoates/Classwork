%% Joshua Oates - Exam 2 - Question 8 9 10
clear all
close all
clc

k = .351;
vimpact = 10;
size = logspace(-6,3,1000);
% size = linspace(1e-6,1e3,600);
rhopart = 3.5;

crater = zeros(1,length(size));
thick = zeros(1,length(size));

for i = 1:length(size)
    crater(i) = dcrater(k,size(i),rhopart,vimpact);
    thick(i) = tcrit(k,size(i),rhopart,vimpact);
end

figure
hold on
plot(size,crater)
% gca().set('xscale','log')
% gca().set('yscale','log')
% legend("Crater Dia")
xlabel("Particle Diameter [cm]")
ylabel("Crater Diameter [cm]")
title("Crater Size vs Particle Diameter")

figure
hold on
plot(size,thick)
gca().set('xscale','log')
gca().set('yscale','log')
% legend("Crater Dia")
xlabel("Particle Mass [g]")
ylabel("Thickness [cm]")
title("Critical Thickness vs Particle Diameter")

disp("The crater diameter of a 1 cm particle at 10 km/s with a density of 3.5 g/cm^3 is: "+string(dcrater(k,1,rhopart,vimpact))+" cm.")
disp("The sheild thickness for a 1 g particle at 10 km/s with a density of 3.5 g/cm^3 is: "+string(tcrit(k,1,rhopart,vimpact))+" cm.")
disp("The plots make sense in so far as the diameter of a crater caused by a particle as well as the thickness of monolithic shielding needed to prevent penetration increases with the size of the particle. Looking at the euqation for crater size, it is sensible that the relation is almost linear, seeing as how the crater diameter is proportional to the particle size raised to the power of 1.056, in other words, almost directly proportional to particle size. The plot of shielding thickness also makes sense because the larger a particle is, the more shielding will be needed.")
%% functions


function out = dcrater(k,dpart,rhopart,vimpact)
    out = k*dpart^1.056*rhopart^.579*vimpact^(2/3);
end

function out = tcrit(k,mpart,rhopart,vimpact)
    out = k*rhopart^(1/6)*mpart^.352*vimpact^.875;
end






