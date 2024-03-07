function [mean] = JoshAvg(input)
%UNTITLED2 Summary of this function goes here
%   takes a 1 dimensional array of numerics and computes the average of it
x = 0;
for i= 1:1:length(input)
    x = x + input(i);
end
mean = x/length(input);
clear x i;
end