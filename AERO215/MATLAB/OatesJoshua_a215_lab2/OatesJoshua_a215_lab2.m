%Joshua Oates aero215 Fall 2021
%Lab2 - flow control
clear
clc
clf
%Part1 - For loop solve dynamic pressure
% q=0.5*density*velocity^2 
% from zero to 250 in steps of 10 m/s at sea level (density = 1.225 kg/m^3)

density_kg_m3 = 1.225;
q_mat = []; %aerodynamic pressure 
v_mat = []; %array version of velocity for plotting

for velocity = 0:10:250 %velocity in m/(s*s) 
    q = 0.5*density_kg_m3*velocity^2;
    v_mat = [v_mat,velocity];
    q_mat = [q_mat,q];
    disp ("v: " + velocity + "m/s^2");
    disp("q: " + q + "N/m^2");
    disp(" ");
end
plot (v_mat,q_mat);


%Part2 - For and if statements

%write for loop that runs 20times, wright if that uses rand() to generate a
%number. print number and say 
%"too big" for greater than or equal to .66
%"just right" for  greater than or equal to 0.33 and less than .66
%"too small" for less than .33

for i = 0:19 %repeat 20 times
    randomNumber = rand(); %create a random number between 0 & 1
    if (randomNumber < .33); %if the number is too small print and end 
        disp("Run#:" + i + "  " + randomNumber + " is too small.");
    elseif (randomNumber < .66); %if the number is not too small or too big print and end
        disp("Run#:" + i + "  " + randomNumber + " is just right!");
    else %if the number is not too small and not correct pring and end
        disp("Run#:" + i + "  " + randomNumber + " is too big.");
    end 
end


