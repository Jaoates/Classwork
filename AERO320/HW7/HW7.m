%% clean up
clear all
close all
clc

%% Solve ODE
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);


tspan = [0,25];
X0 = [1;0];
w=1;

runs = {0;.25;1;2};
n = length(runs);
runs = cat(2,cell(4,1),runs);
runs = cat(2,cell(4,1),runs);
runs = cat(2,runs,{"undamped";"underdamped";"critically damped";"overdamped"});

for i = 1:n
    z = runs{i,3};
    [t,X]=ode45(@eom,tspan,X0,options,z,w);
    runs{i,1} = X;
    runs{i,2} = t;
end
%% plot each

for i = 1:n
    figure
    axe1 = axes;
    title("State Phase "+string(runs{i,4}))
    xlabel("x [m]")
    ylabel("xdot [m/s]")
    axis('equal')
    ylim([-1,1])
    xlim([-1,1])
    hold on


    figure
    axe2 = axes;
    title("x-Time All "+string(runs{i,4}))
    xlabel("time [s]")
    ylabel("x [m]")
    ylim([-1,1])
    hold on

    figure
    axe3 = axes;
    title("xdot-Time "+string(runs{i,4}))
    xlabel("time [s]")
    ylabel("xdot [m/s]")
    ylim([-1,1])
    hold on

    plot(axe1,runs{i,1}(:,1),runs{i,1}(:,2))
    plot(axe2,runs{i,2},runs{i,1}(:,1))
    plot(axe3,runs{i,2},runs{i,1}(:,2))
%     title(axe1,runs{i,4})
%     title(axe2,runs{i,4})
%     title(axe3,runs{i,4})
end



%% plot all 4
figure
axe1 = axes;
title("State Phase All 4")
xlabel("x [m]")
ylabel("xdot [m/s]")
axis('equal')
ylim([-1,1])
xlim([-1,1])
hold on


figure
axe2 = axes;
title("x-Time All 4")
xlabel("time [s]")
ylabel("x [m]")
ylim([-1,1])
hold on

figure
axe3 = axes;
title("xdot-Time All 4")
xlabel("time [s]")
ylabel("xdot [m/s]")
ylim([-1,1])
hold on

for i = 1:n
    plot(axe1,runs{i,1}(:,1),runs{i,1}(:,2))
    plot(axe2,runs{i,2},runs{i,1}(:,1))
    plot(axe3,runs{i,2},runs{i,1}(:,2))
end

legend(axe1,runs{:,4})
legend(axe2,runs{:,4})
legend(axe3,runs{:,4})





%% functions
function Xdot = eom(t,X,z,w)
% d/m = 2wz
% sqrt(k/m) = w
Xdot1 = X(2);
%     Xdot2 = -(d/m)*X(2)-(k/m)*X(1);
Xdot2 = -(2*w*z)*X(2)-(w^2)*X(1);
Xdot = [Xdot1;Xdot2];
end