
clear all
close all
clc



%% [AERO452 PROJECT 1: GOES-18] Josh Oates, David Kumar

load('TLE.mat')

mu = 398600; 
n = tle.MeanMotion;
n = deg2rad(n);
T = 2*pi/n; % sec
MeanAnomaly = tle.MeanAnomaly;

a = 42164000; % km
ecc = 10^-9; % zero
incl = 10^-9; % zero
RAAN = tle.RightAscensionOfAscendingNode;
argp = tle.ArgumentOfPeriapsis;
nu = 0;


[r_vect, v_vect] = keplerian2ijk(a, ecc, incl, RAAN, argp, nu);

r_vect = r_vect/1000;
v_vect = v_vect/1000;

state = [r_vect,v_vect];
t = [0 T];
muearth = 398600;


% figure
% plot3(state(:,1),state(:,2),state(:,3))
% xlabel('x [km]')
% ylabel('y [km]')
% ylabel('z [km]')



%%

addpath("C:\joshFunctionsMatlab\")

% initial relative position and velo/
rinitial = [0 100 0]';
vinitial = [0 0 0]';
r = rinitial;
v = vinitial;

Te = [];
DV = 0;

%% super figures
figure(100)
figure(200)
figure(300)

%% hop one

% xp_2 = 
targetDistance = 40;
dy0 = (norm(r)-targetDistance)/(3*T); % calc hop dv

dv = [0 dy0 0]';
DV = DV+norm(dv);
v=v+dv; % exit hold


t = 0:1:T; % plot vecs
R = zeros([length(t),3]);
V = R;
for i = 1:length(t)
    [R(i,:),V(i,:)]=cwEquations(r,v,t(i),n);
end


[r,v]=cwEquations(r,v,T,n); % complete hop and update locations

dv = -v ;% return to hold
DV = DV+norm(dv);
v = v+dv;
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,state]=ode45(@pig,t,state,options,muearth);


for i =1:length(t)
    dcm = calcDcm(state(i,:));
    Reci(i,:) = dcm*R(i,:)';
end

figure(1)
hold on
% plot3(state(:,1),state(:,2),state(:,3))
plot(R(:,2),R(:,1))
xlabel('vbar [km]')
ylabel('rbar [km]')
% ylabel('z [km]')
title("LVLH Hop 100km - 40km")
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])

figure(2)
hold on
plot(Reci(:,1),Reci(:,2)) % plot hop
% plot3(Reci(:,1),Reci(:,2),Reci(:,3)) % plot hop
xlabel('x [km]')
ylabel('y [km]')
title("rho (ECI) Hop 100km - 40km")
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])

figure(100)
hold on 
plot((t+sum(Te))/(3600*24),vecnorm(R,2,2))
figure(200)
hold on 
plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))
Te = [Te,t(end)]



%% football

a = 40;
dx0 = (a/2)*n;
dv = [1 0 0]'*dx0;
DV = DV+norm(dv);

v = v+dv; % enter football

t = 0:1:T; % plot vecs
R = zeros([length(t),3]);
V = R;
for i = 1:length(t)
    [R(i,:),V(i,:)]=cwEquations(r,v,t(i),n);
end
% plot3(R(:,1),R(:,2),R(:,3)) % plot hop

[r,v]=cwEquations(r,v,T,n); % complete hop and update locations

dv = -v ;% return to hold
DV = DV+norm(dv);
v = v+dv;
t = 0:1:T;


state = state(end,:)';

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,state]=ode45(@pig,t,state,options,muearth);

Reci = [];
for i =1:length(t)
    dcm = calcDcm(state(i,:));
    Reci(i,:) = dcm*R(i,:)';
end

figure(4)
hold on
% plot3(state(:,1),state(:,2),state(:,3))
plot(R(:,2),R(:,1))
xlabel('vbar [km]')
ylabel('rbar [km]')
% ylabel('z [km]')
title("LVLH Football")
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])


figure(5)
hold on
plot(Reci(:,1),Reci(:,2)) % plot hold
axis("equal")
title("rho (ECI) Football")
xlabel('x [km]')
ylabel('y [km]')
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])

figure(100)
hold on 
plot((t+sum(Te))/(3600*24),vecnorm(R,2,2))
figure(200)
hold on 
plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))
Te = [Te,t(end)]


%% hop 2

targetDistance = 1;
dy0 = (norm(r)-targetDistance)/(3*T); % calc hop dv

dv = [0 dy0 0]';
DV = DV+norm(dv);
v=v+dv ;% exit hold

t = 0:1:T; % plot vecs
R = zeros([length(t),3]);
V = R;
for i = 1:length(t)
    [R(i,:),V(i,:)]=cwEquations(r,v,t(i),n);
end


[r,v]=cwEquations(r,v,T,n); % complete hop and update locations

dv = -v; % return to hold
DV = DV+norm(dv);
v = v+dv;
t = 0:1:T;


state = state(end,:)';

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,state]=ode45(@pig,t,state,options,muearth);

Reci = [];
for i =1:length(t)
    dcm = calcDcm(state(i,:));
    Reci(i,:) = dcm*R(i,:)';
end

figure(6)
hold on
% plot3(state(:,1),state(:,2),state(:,3))
plot(R(:,2),R(:,1))
xlabel('vbar [km]')
ylabel('rbar [km]')
% ylabel('z [km]')
title("LVLH Hop 40km - 1km")
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])


figure(7)
hold on
plot(Reci(:,1),Reci(:,2)) % plot hold
axis("equal")
title("rho (ECI) Hop 40km - 1km")
xlabel('x [km]')
ylabel('y [km]')
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])

figure(100)
hold on 
plot((t+sum(Te))/(3600*24),vecnorm(R,2,2))
figure(200)
hold on 
plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))
Te = [Te,t(end)]


%% hold 1
t = 0:1:T;
state = state(end,:)';

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,state]=ode45(@pig,t,state,options,muearth);

% figure
% plot3(state(:,1),state(:,2),state(:,3))
% xlabel('x [km]')
% ylabel('y [km]')
% ylabel('z [km]')
% axis("equal")


Reci = [];
for i =1:length(t)
    dcm = calcDcm(state(i,:));
    Reci(i,:) = dcm*r;
end

figure(8)
hold on
plot(Reci(:,1),Reci(:,2)) % plot hold
axis("equal")
title("rho (ECI) Hold 1km")
xlabel('x [km]')
ylabel('y [km]')
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])

figure(100)
hold on 
plot((sum(Te)+[t(1) t(end)])./(24*3600),[1 1]*norm(r))
figure(200)
hold on 
plot3(r(1)*[1 1],r(2)*[1 1],(sum(Te)+[t(1) t(end)])./(24*3600))
Te = [Te,t(end)]



%% hop 3

targetDistance = .3;
dy0 = (norm(r)-targetDistance)/(3*T); % calc hop dv

dv = [0 dy0 0]';
DV = DV+norm(dv);
v=v+dv; % exit hold
t = 0:1:T; % plot vecs


R = zeros([length(t),3]);
V = R;
for i = 1:length(t)
    [R(i,:),V(i,:)]=cwEquations(r,v,t(i),n);
end



[r,v]=cwEquations(r,v,T,n); % complete hop and update locations

dv = -v; % return to hold
DV = DV+norm(dv);
v = v+dv;

state = state(end,:)';

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,state]=ode45(@pig,t,state,options,muearth);

Reci = [];
for i =1:length(t)
    dcm = calcDcm(state(i,:));
    Reci(i,:) = dcm*R(i,:)';
end

figure(9)
hold on
% plot3(state(:,1),state(:,2),state(:,3))
plot(R(:,2),R(:,1))
xlabel('vbar [km]')
ylabel('rbar [km]')
% ylabel('z [km]')
title("LVLH Hop 1km - 300m")
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])


figure(10)
hold on
plot(Reci(:,1),Reci(:,2)) % plot hold
axis("equal")
title("rho (ECI) Hop 1km - 300m")
xlabel('x [km]')
ylabel('y [km]')
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])

figure(100)
hold on 
plot((t+sum(Te))/(3600*24),vecnorm(R,2,2))
% figure(200)
% hold on 
% plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))
figure(300)
hold on 
plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))
Te = [Te,t(end)]


%% hold 2
t = 0:1:((T/2)+2.105632833061594e+05);
state = state(end,:)';

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,state]=ode45(@pig,t,state,options,muearth);

Reci = [];
for i =1:length(t)
    dcm = calcDcm(state(i,:));
    Reci(i,:) = dcm*r;
end

figure(11)
hold on
plot(Reci(:,1),Reci(:,2)) % plot hold
axis("equal")
title("rho (ECI) Hold 300m")
xlabel('x [km]')
ylabel('y [km]')
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])

figure(100)
hold on 
plot((sum(Te)+[t(1) t(end)])./(24*3600),[1 1]*norm(r))
% figure(200)
% hold on 
% plot3(r(1)*[1 1],r(2)*[1 1],(sum(Te)+[t(1) t(end)])./(24*3600))
figure(300)
hold on 
plot3(r(1)*[1 1],r(2)*[1 1],(sum(Te)+[t(1) t(end)])./(24*3600))

Te = [Te,t(end)]


%% hop 4

targetDistance = .02;
dy0 = (norm(r)-targetDistance)/(3*T); % calc hop dv

dv = [0 dy0 0]';
DV = DV+norm(dv);
v=v+dv; % exit hold
t = 0:1:T; % plot vecs


R = zeros([length(t),3]);
V = R;
for i = 1:length(t)
    [R(i,:),V(i,:)]=cwEquations(r,v,t(i),n);
end



[r,v]=cwEquations(r,v,T,n); % complete hop and update locations

dv = -v; % return to hold
DV = DV+norm(dv);
v = v+dv;

state = state(end,:)';

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,state]=ode45(@pig,t,state,options,muearth);

Reci = [];
for i =1:length(t)
    dcm = calcDcm(state(i,:));
    Reci(i,:) = dcm*R(i,:)';
end

figure(12)
hold on
% plot3(state(:,1),state(:,2),state(:,3))
plot(R(:,2),R(:,1))
xlabel('vbar [km]')
ylabel('rbar [km]')
% ylabel('z [km]')
title("LVLH Hop 300m - 20m")
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])


figure(13)
hold on
plot(Reci(:,1),Reci(:,2)) % plot hold
axis("equal")
title("rho (ECI) Hop 300m - 20m")
xlabel('x [km]')
ylabel('y [km]')
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])

figure(100)
hold on 
plot((t+sum(Te))/(3600*24),vecnorm(R,2,2))
% figure(200)
% hold on 
% plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))
figure(300)
hold on 
plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))
Te = [Te,t(end)]


%% hop error ish thing

targetDistance = .0205;
dy0 = (norm(r)-targetDistance)/(3*T); % calc hop dv

dv = [0 dy0 0]';
DV = DV+norm(dv);
v=v+dv; % exit hold

t = 0:1:T; % plot vecs
R = zeros([length(t),3]);
V = R;
for i = 1:length(t)
    [R(i,:),V(i,:)]=cwEquations(r,v,t(i),n);
end



[r,v]=cwEquations(r,v,T,n); % complete hop and update locations

dv = -v ;% return to hold
DV = DV+norm(dv);
v = v+dv;

state = state(end,:)';

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,state]=ode45(@pig,t,state,options,muearth);

Reci = [];
for i =1:length(t)
    dcm = calcDcm(state(i,:));
    Reci(i,:) = dcm*R(i,:)';
end

figure(14)
hold on
% plot3(state(:,1),state(:,2),state(:,3))
plot(R(:,2),R(:,1))
xlabel('vbar [km]')
ylabel('rbar [km]')
% ylabel('z [km]')
title("LVLH hop error")
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])


figure(15)
hold on
plot(Reci(:,1),Reci(:,2)) % plot hold
axis("equal")
title("rho (ECI) hop error")
xlabel('x [km]')
ylabel('y [km]')
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])

figure(100)
hold on 
plot((t+sum(Te))/(3600*24),vecnorm(R,2,2))

% figure(200)
% hold on 
% plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))
figure(300)
hold on 
plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))

Te = [Te,t(end)]


%% Vbar Approach

dist = norm(r);
timeToDock = 2*3600;
vc = [0,1,0]*(dist/timeToDock);
aCont = [1;0;0] *-2*n*vc;
dvCont = aCont*timeToDock;
dv = norm(vc)+norm(dvCont);
DV = DV+norm(dv);

t = 0:1:timeToDock; % plot vecs
R = ([0;1;0]*linspace(dist,0,length(t)))';


state = state(end,:)';

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,state]=ode45(@pig,t,state,options,muearth);

Reci = [];
for i =1:length(t)
    dcm = calcDcm(state(i,:));
    Reci(i,:) = dcm*R(i,:)';
end

figure(16)
hold on
% plot3(state(:,1),state(:,2),state(:,3))
plot(R(:,2),R(:,1))
xlabel('vbar [km]')
ylabel('rbar [km]')
% ylabel('z [km]')
title("LVLH Approach")
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])


figure(17)
hold on
plot(Reci(:,1),Reci(:,2)) % plot hold
axis("equal")
title("rho (ECI) Approach")
xlabel('x [km]')
ylabel('y [km]')
axis("equal")
scatter(0,0)
legend(["Chaser Trajectory","Target S/C"])
Te = [Te,t(end)]

figure(100)
hold on 
% plot(t+sum(Te),vecnorm(R,2,2))
% 
% figure(200)
% hold on 
% plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))
figure(300)
hold on 
plot3(R(:,1),R(:,2),(t+sum(Te))/(3600*24))

%% other calcs
disp("The time is: "+string(sum(Te))+ " s")
disp("The time is: "+string(sum(Te)/(24*3600))+ " days")
disp("Time needed is: "+string(sum(Te)-(24*3600*10))+"seconds")

syms theta
dist = [100,40,40,1,1,.3,.3,.02,.02,0];
eqn = dist == norm(r_vect)*sin(theta);
sols = [solve(eqn(1)),solve(eqn(2)),solve(eqn(3)),solve(eqn(4)),solve(eqn(5)),solve(eqn(6)),solve(eqn(7)),solve(eqn(8)),solve(eqn(9)),[0;0]];
dTheta = rad2deg(sols(1,:));


thetaA = rad2deg(n*cumsum(Te));
thetaA = [0,thetaA];
thetaB = thetaA+dTheta;

thetaA = mod(thetaA,360);
thetaA = double(thetaA);
thetaB = mod(thetaB,360);
thetaB = double(thetaB);

disp(thetaA')
disp(thetaB')

figure(100)
set(gca,"YScale","log")
title("Distance between Chaser and Target")
xlabel("Time [days]")
ylabel("Distance [km]")

legend(["Hop 100km - 40km","Football","Hop 40km - 1km","Hold 1km","Hop 1km - 300m","Hold 300m","Hop 300m - 20m","Hold 20m"],"Location","best")



% figure(200)
% title("LVLH Whole Trajectory over Time")
% 
% legend(["Hop 100km - 40km","Football","Hop 40km - 1km","Hold 1km","Hop 1km - 300m","Hold 300m","Hop 300m - 20m","Hold 20m","Final Approach"],"Location","best")
% set(gca,"XLim",gca().YLim)

figure(300)
title("LVLH Last Moves over Time")

legend(["Hop 1km - 300m","Hold 300m","Hop 300m - 20m","Hold 20m","Final Approach"],"Location","best")
set(gca,"XLim",gca().YLim)

%% Functions

function [r,dr]=cwEquations(r0,dr0,t,n)
    Phirr = [
        4-3*cos(n*t)    0 0
        6*(sin(n*t)-n*t) 1 0
        0               0 cos(n*t)
    ];

    Phirv = [
        (1/n)*sin(n*t)      (2/n)*(1-cos(n*t))          0
        (2/n)*(cos(n*t)-1)  (1/n)*(4*sin(n*t)-3*n*t)    0
        0                   0                           (1/n)*sin(n*t)
    ];

    Phivr = [
        3*n*sin(n*t)        0   0
        6*n*(cos(n*t)-1)    0   0
        0                   0   -n*sin(n*t)
    ];
    
    Phivv = [
        cos(n*t)    2*sin(n*t)      0
        -2*sin(n*t) 4*cos(n*t)-3    0
        0           0               cos(n*t)
    ];

    r = Phirr * r0 + Phirv * dr0;
    dr= Phivr * r0 + Phivv * dr0;
end

function CE_L = calcDcm(state)
    reci = [state(:,1),state(:,2),state(:,3)]';
    veci = [state(:,4),state(:,5),state(:,6)]';
    x1 = reci/norm(reci);
    h = cross(reci,veci);
    z1 = h/norm(h);
    y1 = cross(z1,x1);
    x2 = [1;0;0];
    y2 = [0;1;0];
    z2 = [0;0;1];
    F1 = [x1,y1,z1];
    F2 = [x2,y2,z2];
    CE_L = joshVectrixDot(F1,F2);
end



