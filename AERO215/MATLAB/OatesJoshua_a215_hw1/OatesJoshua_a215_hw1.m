%Joshua Oates - a125 - Fall 2021
%HW1 - Standard ATM

clear all; close all; clc;

%%part4
%call function and check hand calcs

%set all four hand calc altitudes
% 
% 7,500 ft
% 33,000 ft
% 42,000 ft
% 80,000 ft
testAlt1 = 2286;
testAlt2 = 10058.4;
testAlt3 = 12801.6;
testAlt4 = 24384;

%call function at each hand calc altitude and save them to different
%variables

[TAlt1,PAlt1,rhoAlt1]=stdAtmOatesJoshua(testAlt1);
[TAlt2,PAlt2,rhoAlt2]=stdAtmOatesJoshua(testAlt2);
[TAlt3,PAlt3,rhoAlt3]=stdAtmOatesJoshua(testAlt3);
[TAlt4,PAlt4,rhoAlt4]=stdAtmOatesJoshua(testAlt4);


%%part5
%run function up to 100km

for i=1:1000 % test iterator for 100 different altitudes
    hM(i)=i*100;
    [TM(i),PM(i),rhoM(i)]=stdAtmOatesJoshua(hM(i));
end
hold;
subplot(1,3,1); %set subplot settings
plot(TM,hM); %plot T
ylabel("altitude (m)"); %label altitude
xlabel("temperature (K)"); %label temperature
subplot(1,3,2);
plot(PM,hM);
xlabel("pressure (Pa)");
subplot(1,3,3);
plot(rhoM,hM);
xlabel("density (kg/m^3)");

