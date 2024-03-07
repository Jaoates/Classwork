% AERO 300 - Lab 4 - Joshua Oates
% sections 1 and 2 will use test code provided in class examples, section 3
% will use my own test cases

%% Section 0 - clean up
clear all
close all
clc

%% Section 1 - test Jacobi against same test cases as used in example

disp(" ")
disp("-------------Section 1-------------")
disp(" ")


%For all cases, set the tolearance to
TOL = 0.5 *10^(-6);


% Case 1: A is strictly diagonally row dominant
A = [4, 2;...
     -1, 2];

b = [3; 2];
x_0 = [1; 1];
[x1, X1, k1] = JoshJacobi(A, b, x_0, TOL);
% Case 2: A is NOT strictly diagonally dominant but Jacobi converges


A = [10, -2, -1;...
      -1,  5,  3;...
       2,  2, -2];

b = [1; 2; 3];
x_0 = [1; 1; 1];
[x2, X2, k2] = JoshJacobi(A, b, x_0, TOL);


% Case 3: A is NOT strictly diagonally dominant and Jacobi does not converge
A = [10, -2, -1;...
      -1, -5,  3;...
       2,  2, -2];

b = [1; 2; 3];
x_0 = [1; 1; 1];
[x3, X3, k3] = JoshJacobi(A, b, x_0, TOL);
[x4, X4, k4] = JoshJacobi(A, b, x_0, TOL,200);

str = "In " + k1 + " iterations x1 was found to be: ";
disp(str)
disp(num2str(x1,'%.7f'))

str = "In " + k2 + " iterations x2 was found to be: ";
disp(str)
disp(num2str(x2,'%.7f'))

str = "In " + k3 + " (defualt) iterations x3 was found to not converge, However;";
disp(str)

str = "In " + k4 + " iterations x3 was found to be: ";
disp(str)
disp(num2str(x4,'%.7f'))

%% Section 2 - test GuassSeidel against same test cases as used in example

disp(" ")
disp("-------------Section 2-------------")
disp(" ")

%For all cases, set the tolearance to
TOL = 0.5 *10^(-6);


% Case 1: A is strictly diagonally row dominant
A = [4, 2;...
     -1, 2];

b = [3; 2];
x_0 = [1; 1];
[x1, X1, k1] = JoshGuassSeidel(A, b, x_0, TOL);
% Case 2: A is NOT strictly diagonally dominant but Jacobi converges


A = [10, -2, -1;...
      -1,  5,  3;...
       2,  2, -2];

b = [1; 2; 3];
x_0 = [1; 1; 1];
[x2, X2, k2] = JoshGuassSeidel(A, b, x_0, TOL);


% Case 3: A is NOT strictly diagonally dominant and Jacobi does not converge
A = [10, -2, -1;...
     -1, -5,  3;...
      2,  2, -2];

b = [1; 2; 3];
x_0 = [1; 1; 1];
[x3, X3, k3] = JoshGuassSeidel(A, b, x_0, TOL);

str = "In " + k1 + " iterations x1 was found to be: ";
disp(str)
disp(num2str(x1,'%.7f'))

str = "In " + k2 + " iterations x2 was found to be: ";
disp(str)
disp(num2str(x2,'%.7f'))

str = "In " + k3 + " iterations x3 was found to be: ";
disp(str)
disp(num2str(x3,'%.7f'))

%% Section 3 - test of both against orriginal test cases

disp(" ")
disp("-------------Section 3-------------")
disp(" ")

%For all cases, set the tolearance to
TOL = 0.5 *10^(-6);


% Case 1
A = [0, 0;...
     0, 0];

b = [3; 2];
x_0 = [1; 1];
[x1J, ~, k1J] = JoshJacobi(A, b, x_0, TOL);
[x1G, ~, k1G] = JoshGuassSeidel(A, b, x_0, TOL);


% Case 2
A = [1, 0, 0;...
     0, 1, 0;...
     0, 0, 1];

b = [1; 2; 3];
x_0 = [0; 0; 0];
[x2J, ~, k2J] = JoshJacobi(A, b, x_0, TOL);
[x2G, ~, k2G] = JoshGuassSeidel(A, b, x_0, TOL);


% Case 3
A = [ 10, -2, -1;...
     -1,  3,  3;...
      2,  2, -2];

b = [1; 2; 3];
x_0 = [1; 1; 1];
[x3J, ~, k3J] = JoshJacobi(A, b, x_0, TOL,10000);
[x3G, ~, k3G] = JoshGuassSeidel(A, b, x_0, TOL);

disp(" ")
disp("Case 1")
disp("Jacobi:")
disp("  iterations:")
disp(num2str(k1J))
disp("  x:")
disp(num2str(x1J))
disp(" ")
disp("GuassSeidel:")
disp("  iterations:")
disp(num2str(k1G))
disp("  x:")
disp(num2str(x1G))

disp(" ")
disp("as expected, both failed")

disp(" ")
disp("Case 2")
disp("Jacobi:")
disp("  iterations:")
disp( num2str(k2J))
disp("  x:")
disp( num2str(x2J))
disp(" ")
disp("GuassSeidel:")
disp("  iterations:")
disp(num2str(k2G))
disp("  x:")
disp(num2str(x2G))

disp(" ")
disp("as expected, both quickly succeeded")

disp(" ")
disp("Case 3")
disp("Jacobi:")
disp("  iterations:")
disp( num2str(k3J))
disp("  x:")
disp( num2str(x3J))
disp(" ")
disp("GuassSeidel:")
disp("  iterations:")
disp(num2str(k3G))
disp("  x:")
disp(num2str(x3G))

disp(" ")
disp("although it wasn't garaunteed, both converged. GuassSeidel did so much quicker")

