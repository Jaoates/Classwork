%% Joshua Oates - Exam 2 - Question 2 3 4
clear all
close all
clc

m = logspace(-16,0);

z = 300;
ramOrWake = 1;
flux = zeros(1,length(m));
fluxBack = flux;
for i = 1:length(m)
    flux(i) = Fsp(m(i),z,ramOrWake);
    fluxBack(i) = Fback(m(i));
end

figure
hold on
plot(m,fluxBack)
plot(m,flux)
gca().set('xscale','log')
gca().set('yscale','log')
legend("Background Flux","Actual Flux")
xlabel("Mass [g]")
ylabel("Flux [1/(yr*m^2)]")
title("MM Flux vs Mass")

disp("The addition and reduction factor applied to the background flux is: "+string(fgrav(z)*fshield(z)*fdist(ramOrWake,z))+".")
disp("The flux for mass = 10^-6 g and greater is: "+string(Fsp(10^-6,z,ramOrWake))+" [1/(yr*m^2)].")
disp("This graph makes sense because it is monotomically decreasing with greater masses. This is because the value for flux is for that size mass and larger. In other words as the input mass increases, there are less particles present that are larger than the input mass. The reason that we measure flux of the input mass and all larger masses is because if the input mass would damage the spacecraft, then all larger particles would also cause damage. The adjusted flux being slightly lower makes sense because the spacecraft will be recieving a large amount of shielding from the Earth even though we are studying the ram side of the spacecraft, and there is an increase in density due to Earth's gravity.")

%% functions

function flux = Fsp(m,z,ramOrWake)
flux = Fback(m)*fgrav(z)*fshield(z)*fdist(ramOrWake,z);
end

function out = Fback(m)
out = (3.15576e7)*(A(m)+B(m)+C(m));
end

function out = fgrav(z)
RE = 6378; % km
out = 1 + (RE + 100) / (RE + z);
end

function out = fshield(z)
out = (1+cos(eta(z)))/2;
end

function out = fdist(ramOrWake,z)
% assume ram direction if ramOrWake is 1
arguments
    ramOrWake logical
    z
end
RE = 6378; % km
if ramOrWake
    out = (1.8+3*sqrt(1-((RE+100)/(RE+z))^2))/4;
else
    out = .1;
end
end

function out = eta(z)
RE = 6378; % km
out = asin((RE + 100) / (RE + z));
end

function out = A(m)
out = (15 + (2.2*10^3)*m^(0.306))^(-4.38);
end

function out = B(m)
out = 1.3*10^(-9)*(m + (10^11*m^2) + (10^27*m^4))^(-0.36);
end

function out = C(m)
out = 1.3*10^(-16)*(m + (10^6)*m^2)^(-0.85);
end


