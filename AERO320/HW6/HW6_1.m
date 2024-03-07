% Josh Oates
% HW 6
clear all
J = [...
    [2 -1 0];...
    [-1 3 0];...
    [0 0 1]];

[E,Lam] = eig(J);

disp("Problem 5")
disp("The principle moments of intertia are:")
disp(Lam)
disp("Problem 6")
disp("The components of basis vectors for the principle axes in the refrence fram Fc are:")
disp("X:")
disp(E(:,1))
disp("Y:")
disp(E(:,2))
disp("Z:")
disp(E(:,3))
