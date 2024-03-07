%% housecleaning
clear all
close all
clc
%%single p
X0 = [.5;0];
g = 9.8;
l = 1;
p = .25;
m = 1;
B = .5;
tspan = [0,40];

options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
[t,X]=ode45(@eom,tspan,X0,options,l,p,B,m,g);

% state phase
figure
plot(X(:,1),X(:,2))
xlabel("theta")
ylabel("thetadot")
title("state phase plot, p = .25")
% state time 3d 
% figure
% plot3(X(:,1),X(:,2),t)
% xlabel("theta")
% ylabel("thetadot [1/s]")
% zlabel("time [s]")

% state time 2d
figure
plot(t,X(:,1),t,X(:,2))
legend("theta","thetadot [1/s]")
xlabel("time [s]")
title("state time plot, p = .25")


%%multiple ps
ps = -.25:.125:.25;
n = length(ps);

X0 = [.5;0];
g = 9.8;
l = 1;
m = 1;
B = .5;
tspan = [0,40];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);

%calulate X at each p to test
runs = {};
for i = 1:n
    p = ps(i);
    [t,X]=ode45(@eom,tspan,X0,options,l,p,B,m,g);
    runs{i,1} = X;
    runs{i,2} = t;
end


%%%%%%%%%% phase portrait
figure
hold on
for i = 1:n
%     plot(runs{i,1}(:,1),runs{i,1}(:,2))
    plot3(runs{i,1}(:,1),runs{i,1}(:,2),runs{i,2})
end
xlabel("theta")
ylabel("thetadot [1/s]")
legend("p = "+string(ps))
title("phase portrait, p = var. vals")



%%%%%%%%%% all theta
temp1 = " theta";
temp2 = [];
figure
hold on
for i = 1:n
    plot(runs{i,2},runs{i,1}(:,1));
    temp2 = [temp2,"p = "+string(ps(i))+temp1];
end
legend(temp2)
xlabel("time [s]")
title("theta time plot, p = var. vals")
clear temp2 temp1


%%%%%%%%%% all thetadot
temp1 = " thetadot";
temp2 = [];
figure
hold on
for i = 1:n
    plot(runs{i,2},runs{i,1}(:,2));
    temp2 = [temp2,"p = "+string(ps(i))+temp1];
end
legend(temp2)
xlabel("time [s]")
ylabel("[1/s]")
title("thetadot time plot, p = var. vals")
clear temp2 temp1




%%functions
function Xdot = eom(t,X,l,p,B,m,g);
theta = X(1);
thetadot = X(2);
Xdot = [
    thetadot
    -p*g*(sin(theta))/(l^2+p^2) - B*thetadot/(2*m*(l^2+p^2));
    ];
end