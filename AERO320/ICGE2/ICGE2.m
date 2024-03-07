clear all
close all
clc
addpath("C:\joshFunctionsMatlab\")


%%
% Ir = [[1,0,0];...
%     [0,1,0];...
%     [0,0,0]]*1045.4;
% Ir(3,3) = 19.7;
% 
% Iw = [[1,0,0];...
%     [0,1,0];...
%     [0,0,0]]*.1013;
% Iw(3,3) = .2;

%% no torque
tspan = [0,100];
X0 = zeros(9,1);
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
[t,X] = ode45(@odeFun1,tspan,X0,options);
figure("Name","3")
tiledlayout(2,1)
nexttile
plot(t,X(:,1),t,X(:,2),t,X(:,3))
title("w vs t")
legend("x","y","z")
nexttile
plot(t,X(:,7),t,X(:,8),t,X(:,9))
title("Euler vs t")
legend("phi","theta","psi")

%% torque a
tspan = [0,1000];
X0 = zeros(9,1);
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
[t,X] = ode45(@odeFun2,tspan,X0,options);
figure("Name","4a")
tiledlayout(2,1)
nexttile
plot(t,X(:,1),t,X(:,2),t,X(:,3))
title("w vs t")
legend("x","y","z")
nexttile
plot(t,X(:,7),t,X(:,8),t,X(:,9))
title("Euler vs t")
legend("phi","theta","psi")

%% torque b
tspan = [0,1000];
X0 = zeros(9,1);
%%%%%%%%
X0(3) = .1;
%%%%%%%%
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
[t,X] = ode45(@odeFun3,tspan,X0,options);
figure("Name","4b")
tiledlayout(2,1)
nexttile
plot(t,X(:,1),t,X(:,2),t,X(:,3))
title("w vs t")
legend("x","y","z")
nexttile
plot(t,X(:,7),t,X(:,8),t,X(:,9))
title("Euler vs t")
legend("phi","theta","psi")

%% torque c
tspan = [0,1000];
X0 = zeros(9,1);
%%%%%%%%
% X0(3) = .1;
%%%%%%%%
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
[t,X] = ode45(@odeFun4,tspan,X0,options);
figure("Name","4c")
tiledlayout(2,1)
nexttile
plot(t,X(:,1),t,X(:,2),t,X(:,3))
title("w vs t")
legend("x","y","z")
nexttile
plot(t,X(:,7),t,X(:,8),t,X(:,9))
title("Euler vs t")
legend("phi","theta","psi")

%% torque d
% tspan = [0,1000]
X0 = zeros(9,1);
%%%%%%%%
X0(3) = .1;
%%%%%%%%
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
[t,X] = ode45(@odeFun4,tspan,X0,options);
figure("Name","4d")
tiledlayout(2,1)
nexttile
plot(t,X(:,1),t,X(:,2),t,X(:,3))
title("w vs t")
legend("x","y","z")
nexttile
plot(t,X(:,7),t,X(:,8),t,X(:,9))
title("Euler vs t")
legend("phi","theta","psi")


%% functions

function dX = odeFun1(t,X)
w = X(1:3);
wrel = X(4:6);
E = X(7:9);

%%%%%%%
dwrel =[0;0;.05];
%%%%%%%
Td = [0;0;0];
%%%%%%%

Ir = [[1,0,0];...
    [0,1,0];...
    [0,0,0]]*1045.4;
Ir(3,3) = 19.7;

Iw = [[1,0,0];...
    [0,1,0];...
    [0,0,0]]*.1013;
Iw(3,3) = .2;

mat = [[cos(E(2)),      sin(E(2))*sin(E(1)),            sin(E(2))*cos(E(1))];...
    [0,                 cos(E(2))*cos(E(1)),           -cos(E(2))*sin(E(1))];...
    [0,                 sin(E(1)),                      cos(E(1))           ]]*(1/cos(E(2)));

wx = joshCross(w);
ww=(w+wrel);
wwx = joshCross(ww);
dw = -inv(Ir+Iw)*(wx*Ir*w+Iw*dwrel+wwx*Iw*ww+Td);

Edot = mat*w;

dX = [dw;dwrel;Edot];
end

function dX = odeFun2(t,X)
w = X(1:3);
% wrel = X(4:6);
wrel = zeros(3,1);
E = X(7:9);

%%%%%%%
dwrel =[0;0;0];
%%%%%%%
Td = [.1;0;0];
%%%%%%%

Ir = [[1,0,0];...
    [0,1,0];...
    [0,0,0]]*1045.4;
Ir(3,3) = 19.7;

Iw = [[1,0,0];...
    [0,1,0];...
    [0,0,0]]*.1013;
Iw(3,3) = .2;

mat = [[cos(E(2)),      sin(E(2))*sin(E(1)),            sin(E(2))*cos(E(1))];...
    [0,                 cos(E(2))*cos(E(1)),           -cos(E(2))*sin(E(1))];...
    [0,                 sin(E(1)),                      cos(E(1))           ]]*(1/cos(E(2)));

wx = joshCross(w);
ww=(w+wrel);
wwx = joshCross(ww);
dw = inv(Ir+Iw)*(wx*Ir*w+Iw*dwrel+wwx*Iw*ww+Td);

Edot = mat*w;

dX = [dw;dwrel;Edot];
end

function dX = odeFun3(t,X)
w = X(1:3);
% wrel = X(4:6);
wrel = zeros(3,1);
E = X(7:9);

%%%%%%%
dwrel =[0;0;0];
%%%%%%%
Td = [.1;0;0];
%%%%%%%

Ir = [[1,0,0];...
    [0,1,0];...
    [0,0,0]]*1045.4;
Ir(3,3) = 19.7;

Iw = [[1,0,0];...
    [0,1,0];...
    [0,0,0]]*.1013;
Iw(3,3) = .2;

mat = [[cos(E(2)),      sin(E(2))*sin(E(1)),            sin(E(2))*cos(E(1))];...
    [0,                 cos(E(2))*cos(E(1)),           -cos(E(2))*sin(E(1))];...
    [0,                 sin(E(1)),                      cos(E(1))           ]]*(1/cos(E(2)));

wx = joshCross(w);
ww=(w+wrel);
wwx = joshCross(ww);
dw = -inv(Ir+Iw)*(wx*Ir*w+Iw*dwrel+wwx*Iw*ww+Td);

Edot = mat*w;

dX = [dw;dwrel;Edot];
end

function dX = odeFun4(t,X)
w = X(1:3);
% wrel = X(4:6);
wrel = zeros(3,1);
wrel(3) = 100;
E = X(7:9);

%%%%%%%
dwrel =[0;0;0];
%%%%%%%
Td = [.1;0;0];
%%%%%%%

Ir = [[1,0,0];...
    [0,1,0];...
    [0,0,0]]*1045.4;
Ir(3,3) = 19.7;

Iw = [[1,0,0];...
    [0,1,0];...
    [0,0,0]]*.1013;
Iw(3,3) = .2;

mat = [[cos(E(2)),      sin(E(2))*sin(E(1)),            sin(E(2))*cos(E(1))];...
    [0,                 cos(E(2))*cos(E(1)),           -cos(E(2))*sin(E(1))];...
    [0,                 sin(E(1)),                      cos(E(1))           ]]*(1/cos(E(2)));

wx = joshCross(w);
ww=(w+wrel);
wwx = joshCross(ww);
dw = -inv(Ir+Iw)*(wx*Ir*w+Iw*dwrel+wwx*Iw*ww+Td);

Edot = mat*w;

dX = [dw;dwrel;Edot];
end