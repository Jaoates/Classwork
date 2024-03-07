%% Question 1
clear all
close all
clc

rho = [2300,2700,1700];%kg/m^3
rho = rho*1e3;%g/m^3
rho = rho/(1e2)^3;% g/cm^3
si = [1.567, 1.118e2, 1.024e3]; %MeV cm^2/g
al = [1.518, 1.095e2, 9.859e2];
gr = [1.619, 1.409e2, 1.404e3];

ke = 2;%Mev

Rsiw=ke./si;
Ralw=ke./al;
Rgrw=ke./gr;

Rsi=ke./(rho(1)*si);
Ral=ke./(rho(2)*al);
Rgr=ke./(rho(3)*gr);

S = [si;al;gr];
Rw = [Rsiw;Ralw;Rgrw];
R = [Rsi;Ral;Rgr];

label1 = ["Si","Al","gr"];
label2 = ["Electron","Protons","Alpha"];

disp("Penetration Depth [cm]:")
spec = '\t %s \t\t\t %s \t\t\t %s\n';
fprintf(spec,label2)
spec = '%s \t %d \t\t %d \t\t %d\n';
for i = 1:3
    fprintf(spec,label1(i),R(i,:))
end



disp("Specific Penetration Depth [cm^2/g]:")
spec = '\t %s \t\t\t %s \t\t\t %s\n';
fprintf(spec,label2)
spec = '%s \t %d \t\t %d \t\t %d\n';
for i = 1:3
    fprintf(spec,label1(i),Rw(i,:))
end

disp("Stopping Power:")
spec = '\t %s \t\t\t %s \t\t\t %s\n';
fprintf(spec,label2)
spec = '%s \t %d \t\t %d \t\t %d\n';
for i = 1:3
    fprintf(spec,label1(i),S(i,:))
end

disp("For all cases, electrons are the most likely to do damage and alpha particles are the least likely.")
disp("Since gr has the lowest penetration depth for a given mass, it is the best sheiliding material for space applications where mass is the most important parameter to minimize.")
