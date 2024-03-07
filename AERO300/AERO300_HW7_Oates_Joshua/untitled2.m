close all
clear all
clc
figure
f = @(t) plot(t,1,'r*');
fanimator(f)
axis equal
hold off
playAnimation