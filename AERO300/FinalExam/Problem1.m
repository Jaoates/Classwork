%joshua Oates
%Final Exam
%% zero
clear all;
close all;
clc;
%% one
c = .5;
f = @(x)exp(-100*(x-1).^2);
xdom = [0,2];
tdom = [0,5];
xnum = 200;
tnum = 1000;
[X,T,W1] = advectionEquation1DFTCS(c,f,xdom,tdom,xnum,tnum);


figure
plot(X,W1(:,end))
title("Advection 1D - FTCS - T=5")
xlabel("X")
ylabel("advection")

disp("1.f) It looks unstable, similar to unstable heat graphs.")


%% two
c = .5;
f = @(x)exp(-100*(x-1).^2);
xdom = [0,2];
tdom = [0,5];
xnum = 200;
tnum = 1000;
[X,T,W2] = advectionEquation1DLF(c,f,xdom,tdom,xnum,tnum);

k = 5/1000;
times = [0,.5,1,1.5,5];
times = times./k;

figure
hold on
for i = 1:4 % change this value from 1:4 to 1:5 to see first 4 times vs T=5
    plot(X,W2(:,times(i)+1))
end
title("Advection 1D - LF")
xlabel("X")
ylabel("advection")
legend("T=0","T=.5","T=1","T=1.5","T=5")
disp("2.e) It looks as though even though this scheme is stable, the diffusity issue causes decay in amplitude over time")

%% three
c = .5;
f = @(x)exp(-100*(x-1).^2);
xdom = [0,2];
tdom = [0,5];
xnum = 200;
tnum = 1000;
[X,T,W3] = advectionEquation1DLW(c,f,xdom,tdom,xnum,tnum);

k = 5/1000;
times = [0,.5,1,1.5,5];
times = times./k;

figure
hold on
for i = 1:4 % change this value from 1:4 to 1:5 to see first 4 times vs T=5
    plot(X,W3(:,times(i)+1))
end
title("Advection 1D - LW")
xlabel("X")
ylabel("advection")
legend("T=0","T=.5","T=1","T=1.5","T=5")
disp("3.e) It looks as though the diffusivity issue is solved for most values less than T=5. at T=5 it is much reduced from the previous scheme. This is obvious whe comparing figures 2 and 3 since the amplitude in figure 3 does not noticably change over time until T=5. This is very distinctly different from figure 2.")


