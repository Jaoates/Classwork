%% Hw1 - A402
clear all
close all
clc
%% problem 1



%% problem 2
clear all
close all
mdotp = 650; %kg/s      mass flow propelent
c = 4400; % m/s         effective velocity
mStruct = 21500; %kg    
mProp = 180000; %kg
mPay = 15000; %kg
m = mStruct + mProp + mPay;
g = 9.807;

F = mdotp*c
Isp = F/(mdotp*g)
I = F*(mProp/mdotp)
T_W = F/(m*g)
mInert = mStruct;
mI_m = mInert/m
Dv = -c*log((mInert+mPay)/m)

%% problem 3
clear
close all

g = 9.807;
Isp = [60;320;950;3000];
Dv= @(fInert) Isp.*log(1/fInert).*g
fInert = 0:.001:1;

% DV = []

for i = 1:length(fInert)
    DV(:,i) = Dv(fInert(i));
end

figure
hold on 
% tiledlayout(2,2)

for i = 1:4
    % nexttile
    plot(fInert,DV(i,:))
    set(gca,"yscale","log")
    xlabel("Inert Fraction")
    ylabel("Delta V")
end
legend(string(Isp)+"s")

fInert = [.05,.4,.05,.05];
TargetDV = [50,3260+680,3210+640+2940,(3210+4500+10230)*2];

names= ["a","b","c","d"]';

for i = 1:length(fInert)
    disp("part:"+ string(names(i)))
    disp("the following Isp values would work for this situation:")
    disp(string(Isp(Dv(fInert(i))>=TargetDV(i)))+"s")
    scatter(fInert(i),TargetDV(i))
end

title("Feasibility lines for 4 missions")
legend([[string(Isp)+"s"],names],"location","best")





