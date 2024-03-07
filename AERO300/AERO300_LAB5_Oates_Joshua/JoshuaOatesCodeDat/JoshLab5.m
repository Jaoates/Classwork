% Joshua Oates - Lab 5

%% Section 0 - clean up
clear all;
close all;
clc;
addpath('Data')
%% Section 1 - plot the wing

air1= load('airfoil_1.txt'); % load dat
air2= load('airfoil_2.txt');

p1 = polyfit(air1(:,1),air1(:,2),5); % create best fits
p2 = polyfit(air2(:,1),air2(:,2),5);
air1Chord = abs([min(air1(:,1)),max(air1(:,1))]); % create range that we plot over
air2Chord = abs([min(air2(:,1)),max(air2(:,1))]);
X1 = linspace(air1Chord(1),air1Chord(2)); % create linspace to plot over, X is chord axis
X2 = linspace(air2Chord(1),air2Chord(2));
Z1 = polyval(p1,X1); % create Z vectors to plot over Z is height axis
Z2 = polyval(p2,X2);

% figure
% hold on
% plot(air1(:,1),air1(:,2),'* r')
% plot(air2(:,1),air2(:,2),'* b')
% plot(X1,Z1,'r')
% plot(X2,Z2,'b')
% axis('equal')

Y1 = zeros(length(X1),1); % create Y vectors to plot over Y is wingspan axis
Y2 = ones(length(X1),1);
X = [X1',X2'];
Z = [Z1',Z2'];
Y = [Y1,Y2];
figure
hold on 
surf(X,Y,Z,EdgeColor = 'none') % plot + and - Z over domains
surf(X,Y,-Z,EdgeColor = 'none')
axis('equal')
title("3d Airfoil")
xlabel("x (chords)")
ylabel("y (chords)")
zlabel("z (chords)")
campos([-2,2,1.5]) % change camera so that the 3d nature of th wing plot is shown
camtarget([1,0,0])

clear air1 air2 air2Chord air1Chord X1 X2 X Y1 Y2 Y Z1 Z2 Z Ydomain p1 p2

%% Section 2 - plot wind tunnel data
cl = load('cl_data.txt'); % load dat
cd = load('cd_data.txt');
cltrim = [cl((1:13),(1:2));cl((36:49),(1:2))]; % trim out data that will not make a nice line

p1 = polyfit(cltrim(:,1),cltrim(:,2),1); % create line of best fit
X1 = linspace(-6,9); % create plotable data
Y1 = polyval(p1,X1);  

p2 = polyfit(cd(:,1),cd(:,2),2); % create quadratic of best fit
X2 = linspace(min(cd(:,1)),max(cd(:,1))); % create plotable data
Y2 = polyval(p2,X2);


figure
hold on 
plot(cl(:,1),cl(:,2),'. g') 
plot(cltrim(:,1),cltrim(:,2),'. b')
plot(cd(:,1),cd(:,2),'. r')
plot(X2,Y2, 'r')
plot(X1,Y1, 'k')
plot(0,p1(2), '* k')
legend("Coeff Lift", "Cl Used for Best Fit","Coeff Drag","Cd Best Fit","Cl Best Fit","Cl0",Location='northwest')
title("Plot of Aerodynamic Properties vs Angle")
xlabel("Angle of Attack")

disp("I found Cl0 to be "+ p1(2))
disp("I found the proportionality constant between angle of attack and lift coefficient to be "+ p1(1))
disp("I found the equation of the quadratic of best fit for the Coeff Drag to be: "+p2(1)+ "X^2+ "+p2(2)+ "X+ "+p2(3))






