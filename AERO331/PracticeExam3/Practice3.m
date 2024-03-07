%% Practice Exam 2
clear all
close all
clc

% Ei = [30e6 10e6];
Ei = [1 1 1 1];
% alpha = [6.e-6 13e-6];
E1 = 1;
Ei_E1 = Ei./E1;

syms t d

Ai = [[1 1]*t*d,[1 1]*t*d];

zip = [[0 0 1 -1]*d];
yip = [[1 -1 0 0]*d];

% b = 1.5;
% a1 = 1;
% a2 = 2;
% 
% Izoi = [a1^3*b/12 a2^3*b/12];
% Iyoi = [a1*b^3/12 a2*b^3/12];

Izoi = [t^3*d t^3*d d^3*t d^3*t]./12
Iyoi = [d^3*t d^3*t t^3*d t^3*d]./12

Iyzoi = [1 1 1 1]*t*d*(t^2+d^2)/12;

thing = joshAdvBeam(Ai,yip,zip,Iyoi,Izoi,Iyzoi,Ei_E1);









