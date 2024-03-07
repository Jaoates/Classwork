clear all
close all
d = importdata("A356_Lab1_Cylinders_Run1_data.csv");
d = d.data;
d1.t = d(:,1);
d1.T1 = d(:,2);
d1.T2 = d(:,3);
d1.T3 = d(:,4);
d1.T4 = d(:,5);

d = importdata("A356_Lab1_Cylinders_Run2_data.csv");
d = d.data;
d2.t = d(:,1);
d2.T1 = d(:,2);
d2.T2 = d(:,3);
d2.T3 = d(:,4);
d2.T4 = d(:,5);
clear d

pamb1 = polyfit(d1.t,d1.T4,1);
pamb2 = polyfit(d2.t,d2.T4,1);

pamb1 = @(t) pamb1(1)*t+pamb1(2);
pamb2 = @(t) pamb2(1)*t+pamb2(2);

d2.aveWin = [200,350];
d2.aveWin = d2.t>d2.aveWin(1) & d2.t<d2.aveWin(2);
d2.dT1 = polyfit(d2.t(d2.aveWin),d2.T1(d2.aveWin),1);
d2.dT2 = polyfit(d2.t(d2.aveWin),d2.T2(d2.aveWin),1);
d2.dT3 = polyfit(d2.t(d2.aveWin),d2.T3(d2.aveWin),1);
d2.dT4 = polyfit(d2.t(d2.aveWin),d2.T4(d2.aveWin),1);

d1.aveWin = [400,500];
d1.aveWin = d1.t>d1.aveWin(1) & d1.t<d1.aveWin(2);
d1.dT1 = polyfit(d1.t(d1.aveWin),d1.T1(d1.aveWin),1);
d1.dT2 = polyfit(d1.t(d1.aveWin),d1.T2(d1.aveWin),1);
d1.dT3 = polyfit(d1.t(d1.aveWin),d1.T3(d1.aveWin),1);
d1.dT4 = polyfit(d1.t(d1.aveWin),d1.T4(d1.aveWin),1);


i = 1;
c(i).material = "Al";
c(i).color = "White";
c(i).nomD = ".75";
c(i).m = 29.431;
c(i).d = .7556;%in
c(i).L = 1.509;
c(i).c = 0.89;% J/g
c(i).dT = d1.dT1(1);

i = 2;
c(i).material = "Al";
c(i).color = "Polished";
c(i).nomD = ".75";
c(i).m = 27.743;
c(i).d = .7335;
c(i).L = 1.50;
c(i).c = 0.89;% J/g
c(i).dT = d1.dT2(1);


i = 3;
c(i).material = "Al";
c(i).color = "Black";
c(i).nomD = ".75";
c(i).m = 29.518;
c(i).d = .757;
c(i).L = 1.507;
c(i).c = 0.89;% J/g
c(i).dT = mean([d1.dT3(1),d2.dT3(1)]);


i = 4;
c(i).material = "Al";
c(i).color = "Black";
c(i).nomD = "1.0";
c(i).m = 52.705;
c(i).d = 1.008;
c(i).L = 1.512;
c(i).c = 0.89;% J/gK
c(i).dT = d2.dT1(1);


i = 5;
c(i).material = "Br";
c(i).color = "Black";
c(i).nomD = ".75";
c(i).m = 91.518;
c(i).d = .7545;
c(i).L = 1.511;
c(i).c = 0.380;% J/gK
c(i).dT = d2.dT2(1);

clear i

for i = 1:length(c)
c(i).d = c(i).d*0.0254;%in -> m
c(i).L = c(i).L*0.0254;
end

%%%%%% 
% AL W .75
% AL S .75
% AL B .75
leg = c(1).material +" " +c(1).color +" " + c(1).nomD;
leg = [leg,(c(2).material +" " +c(2).color +" " + c(2).nomD)];
leg = [leg,(c(3).material +" " +c(3).color +" " + c(3).nomD)];

leg = [leg,"Ambient","Ambient Best Fit"];
f2 = figure 
hold on
plot(d1.t,d1.T1,d1.t,d1.T2,d1.t,d1.T3,d1.t,d1.T4)
fplot(pamb1)
legend(leg,Location="best")
title("Run 1")
ylabel("Temp [K]")
xlabel("Time [s]")

% AL B 1
% BR S 1
% AL B .75
leg = c(4).material +" " +c(4).color +" " + c(4).nomD;
leg = [leg,(c(5).material +" " +c(5).color +" " + c(5).nomD)];
leg = [leg,(c(3).material +" " +c(3).color +" " + c(3).nomD)];
leg = [leg,"Ambient","Ambient Best Fit"];

f1 = figure
hold on
plot(d2.t,d2.T1,d2.t,d2.T2,d2.t,d2.T3,d2.t,d2.T4)
fplot(pamb2)

legend(leg,Location="best")
title("Run 2")
ylabel("Temp [K]")
xlabel("Time [s]")

% find the Solar constant


I = @(c)(c.m*c.c*c.dT)/(c.d*c.L);

for i = 1:length(c)
    c(i).I = I(c(i));
end

exportgraphics(f1,"CylindersRun1.png","Resolution",300)
exportgraphics(f2,"CylindersRun2.png","Resolution",300)

