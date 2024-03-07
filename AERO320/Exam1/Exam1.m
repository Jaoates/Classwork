%% Midterm exam 1- Joshua Oates
clear all
close all
clc

addpath('C:\joshFunctionsMatlab\')


%% Problem 2

[Cx,Cy,Cz] = joshAxisRotation();

syms a t

c21=Cz(a)
