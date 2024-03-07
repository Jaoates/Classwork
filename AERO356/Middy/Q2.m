%% Question 2
clear all
close all
clc

label = ["photon","electron","proton"];
abDos0 = [6,5,5]*1e-3;
abDos1 = [18,15,8]*1e-3;
abDos2 = [21,19,9]*1e-3; %  /m

scaleFactor = 1;%hr
yr = 365.25*scaleFactor; %hr
m = yr/12;
hrsInMonth = m*scaleFactor
abDos0 = abDos0/m;
abDos1 = abDos1/m;
abDos2 = abDos2/m;

n = round(yr*3);
D = ones(n,1)*abDos0;
D(round(4*m):round(5*m)-1,:) =  ones(round(5*m)-round(4*m),1)*abDos1;
D(round(11*m):round(12*m)-1,:) =  ones(round(12*m)-round(11*m),1)*abDos2;
dD = 3*abDos0./(6*m);
for i = 22*m:28*m
    D(round(i),:) = D(round(i-1),:)+dD;
end
cum = cumsum(D);
Wr = [1,1,2];
cumE = Wr.*cum;

abDos = sum(cum,2);
abDosE = sum(cumE,2);


for i = 1:3*12
    mDos(i) = sum(sum(D(round(m*(i-1)+1):round(m*(i+0)),:),2));
    mDosE(i) = sum(sum(Wr.*D(round(m*(i-1)+1):round(m*(i+0)),:),2));
end

months = linspace(1,3*12,n);

figure
hold on
plot(months,abDos)
% plot(months,abDosE)    
% bar(mDosE)
bar(mDos)
% plot([0,months(end)],[1,1]*.9)
legend("Absorbed","Monthly Absorbed",Location="best")
title("Cumulative 3 year")

figure
hold on
plot(months(1:round(length(abDosE)/3)),abDosE(1:round(length(abDosE)/3)))
plot(months(1:round(length(abDosE)/3)+2),-abDosE(round(length(abDosE)/3)) + abDosE(round(length(abDosE)/3):round(length(abDosE)*2/3)))
plot(months(1:round(length(abDosE)/3)+1),-abDosE(round(length(abDosE)*2/3)) +abDosE(round(length(abDosE)*2/3:end)))
plot([0,months(round(length(abDosE)/3))],[1,1]*.5)
legend("year 1","year 2","year 3","Year Limit",Location="best")
title("Cumulative 3 year Effective vs Absorbed")
ylabel("Gy/Sv")
xlabel("Months")

figure
hold on
plot(months,abDos)
plot(months,abDosE)    
bar(mDosE)
bar(mDos)
plot([0,months(end)],[1,1]*.9)
legend("Absorbed","Effective","Monthly Absorbed","Monthly Effective","Lifetime Limit",Location="best")
title("Cumulative 3 year Effective vs Absorbed")
ylabel("Gy/Sv")
xlabel("Months")
