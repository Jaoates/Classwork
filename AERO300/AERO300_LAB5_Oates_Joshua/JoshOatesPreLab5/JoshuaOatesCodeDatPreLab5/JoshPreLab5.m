% Joshua Oates
% this is my code for prelab 5

%% Section 0 - clean up
clear all;
close all;
clc;

%% Section 1
% NOTE: This must be done before lab. If it is not completed, you will not be allowed into lab. 
% Use the polyfit() command to find the coefficients of a 3rd order polynomial that passes through the 
% points (-5,4), (-2, -1), (4, 2), and (5,-5).  
x = [-5,-2,4,5];
y = [4,-1,2,-5];

p = polyfit(x,y,3);

% Evaluate and plot the polynomial on the interval [-5:.1:5] using
% the polyval() command.  Also, plot the points on the line from above using red x’s as data markers. 

hold on
plot(x,y,'x r')
X = -5:.1:5;
Y = polyval(p,X);
plot(X,Y,'b')
xlabel("x")
ylabel("y")
title("polyfit of given data")
legend("Given Data","Polyfit Degree - 3")
% Be sure to label your plot.  Finally, comment on how is the polyfit() command similar to Newton's Divided 
% Difference method?
disp("polyfit() is similar to Newton's Divided Difference method in that they both aim to interpolate, or get as close as possible, to the given data")
