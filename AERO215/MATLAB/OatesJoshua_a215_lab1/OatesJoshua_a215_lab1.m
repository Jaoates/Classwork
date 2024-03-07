% Joshua Oates Aero215 - 03
% Lab 1 Law of Cosines
% 
%set length of sideA and sideB
sideA = 10;
sideB = 10;
%set angleAB in degrees
angleAB = 60;
%create output var and call function to test the lawCos function
sideC = lawCos(sideA,sideB,angleAB);
%send the output and inputs to command window
clc;
disp("sides:");
disp("A: " + sideA);
disp("B: " + sideB);
disp("C: " + sideC);

