% constants for HW5
clear all
close all
clc

syms y z
E1 = 30e6;
E2 = 10e6;
a1 = 6.5e-6;
a2 = 13e-6;
DT = 10*y^2+y^3;

PT = int(int(E1*a1*DT,y,-1.1,-.1),z,-.75,.75)+int(int(E2*a2*DT,y,-.1,1.9),z,-.75,.75);

MyT = int(int(E1*a1*DT*z,y,-1.1,-.1),z,-.75,.75)+int(int(E2*a2*DT*z,y,-.1,1.9),z,-.75,.75);

MzT = int(int(E1*a1*DT*y,y,-1.1,-.1),z,-.75,.75)+int(int(E2*a2*DT*y,y,-.1,1.9),z,-.75,.75);

PT = double(vpa(PT));
MyT = double(vpa(MyT));
MzT = double(vpa(MzT));

disp("If the temperature is applied then: ")
disp("PT = "+string(PT)+" lbs")
disp("MyT = "+string(MyT)+" lbs-in")
disp("MzT = "+string(MzT)+" lbs-in")