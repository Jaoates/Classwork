clear all
close all
fclose("all")
clc

X = fopen("joshAstroConstantsTable.txt","r");

x = fgetl(X);
lines=[];
while x ~= -1
    lines = [lines;string(x)];
    x = fgetl(X);
end
fclose("all")

lines = char(lines);
[n,m] = size(lines);

lsplits = {};
lens = [];
for i = 1:n
    lsplit = split(lines(i,:));
    lens = [lens;length(lsplit)];
    lsplits{i} = lsplit'
end
lsplits = lsplits'
T = cell2table(lsplits)
clear lsplits lsplit lens lines x X m
%%
clear all
T = load("intermediateTable.mat");
T = T.T;
T = string(T);
T = T(:,1:end-1);
T = table(T);
% S = T('Sun',:)

%%
T2=readtable("joshAstroConstantsTable.txt")

