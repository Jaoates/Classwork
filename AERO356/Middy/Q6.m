%% Question 6
clear all
close all
clc

disp("from Ix = I0*exp(-mu*x)")
disp("we get mu = (-1/x) *log(CPM/CPM0)")

CPM = 1984/30;
CPM0 = 2439/10;
x = 3;
mu = (-1/x) *log(CPM/CPM0)
