clear all
close all
clc

q = 1.60217663e-19;
kb = 5.6704e-8;
mi = 6.6335209e-26;
me = 9.1093837e-31;

d = .02;
l = .075;
d = d*0.0254;
l = l*0.0254;
c = pi*d;
Ap = c*l+pi*(d/2)^2

% Ap = 6.2e-6;

d = readmatrix("LangmuirProfessional.xlsx");

% Plot input data
figure
hold on
plot(d(:,1),d(:,2),".")
xlabel("V")
ylabel("I")
title("I-V curve")

%%%%%%%%%%%%%%%%%%%%%%%% region II
% select and plot range we take to be linear
linParam = [20,40];% volts
% xline(linParam)

% calc best fit
[~,linRange] = min(abs(d(:,1)-linParam)); %index
p1 = polyfit(d(linRange(1):linRange(2),1) , d(linRange(1):linRange(2),2) , 1);

linY = p1(end)+p1(1).*linParam;
Te = (linParam(1)-linParam(2))/(log(linY(1)/linY(2)));
Te = Te*q/kb

% plot(linParam,linY,"*")
% plot the linear section
linParam = [10,50]; % volts
linY = p1(end)+p1(1).*linParam;
plot(linParam,linY)

%%%%%%%%%%%%%%%%%%%%%%%% region I
% select and plot range we take to be linear
linParam = [-30,0];% volts

% calc best fit
[~,linRange] = min(abs(d(:,1)-linParam)); %index
p2 = polyfit(d(linRange(1):linRange(2),1) , d(linRange(1):linRange(2),2) , 1);

linY = p2(end)+p2(1).*linParam;
% Te = (linParam(1)-linParam(2))/(log(linY(1)/linY(2)))

% plot the linear section
linParam = [-35,20]; % volts
linY = p2(end)+p2(1).*linParam;
plot(linParam,linY)

%%%%%%%%%%%%%%%%%%%%%%%% region III
% select and plot range we take to be linear
linParam = [60,70];% volts

% calc best fit
[~,linRange] = min(abs(d(:,1)-linParam)); %index
p3 = polyfit(d(linRange(1):linRange(2),1) , d(linRange(1):linRange(2),2) , 1);

linY = p3(end)+p3(1).*linParam;
% Te = (linParam(1)-linParam(2))/(log(linY(1)/linY(2)))

% plot the linear section
linParam = [35,70]; % volts
linY = p3(end)+p3(1).*linParam;
plot(linParam,linY)

%%%%%%%%%%%%%%%%%%%%%%%% sat vals
syms x
y1 = p1(end)+p1(1).*x;
y2 = p2(end)+p2(1).*x;
y3 = p3(end)+p3(1).*x;

sati(1) = double(solve(y1==y2)); % x
sate(1) = double(solve(y1==y3));
sati(2) = double(subs(y1,sati(1))); % y
sate(2) = double(subs(y3,sate(1)));
plot(sati(1),sati(2),"*")
plot(sate(1),sate(2),"*")

disp("Phi_p is "+string(sate(1))+" V")
[~,If]=min(abs(d(:,2)));
Vf = d(If,1);
disp("Floating voltage (Vf) is "+string(Vf)+" V")

syms ni ne
eqn = -sati(2) == .6*q*ni*sqrt(kb*Te/mi)*Ap;
ni = double(solve(eqn,ni))


Ve = sqrt(8*kb*Te/(pi*me))
eqn = sate(2) == .25*q*ne*Ve*Ap;
ne = double(solve(eqn,ne))


%%%%%%%%%%%%%%%%%%%%%%%% legend and labels

% labels = ["Measurements" "Linear Region Bounds" "Linear Region Bounds" "Region II Linear" "Region I Linear" "Region III Linear" "Ion Saturation" "Electron Saturation"];
labels = ["Measurements" "Region II Linear" "Region I Linear" "Region III Linear" "Ion Saturation" "Electron Saturation"];
legend(labels,"location","best")







