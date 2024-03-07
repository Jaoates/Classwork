%% HW7 - Joshua Oates
clear all
close all
clc

y=load('assignment_7_data.mat');
y=y.data;

%% Plot Raw Data
figure 
hold on
plot(y)


%% Kalman Inputs

F = 1
G = 1
H = 1
Q = eye(1)*0
R = eye(1)*.01
L = []

P0 = 1
x0 = y(1)

P = [P0];
xh = [x0];

u = 0;
u = u*ones(1,length(y));
%% Kalman Filter

for k = 2:length(y)
    xhk_1 = xh(:,:,k-1);
    Pk_1 = P(:,:,k-1);
    [xhk,yhk,Pk] = Kfilter(xhk_1,Pk_1,u(k),y(k),F,G,H,Q,R);
    xh(:,:,k) = xhk;
    yh(:,:,k) = yhk;
    P(:,:,k) = Pk;
end

%% plot kalman outputs
% plot(squeeze(yh))
plot(squeeze(xh))
legend(["Signal","Kalman Filter Estimation"])
xlabel("k")
ylabel("V")
title("Kalman Plot")


%% Functions
function [xhk,yhk,Pk] = Kfilter(xhk_1,Pk_1,u,y,F,G,H,Q,R)
xhmk = F*xhk_1+G*u;
Pkm = F*Pk_1*F'+Q;

Wk = H*Pkm*H'+R;
Kk = Pkm*H'/Wk;
yhk = H*xhmk;
xhk = xhmk+Kk*(y-yhk);
Pk = Pkm-Kk*H*Pkm'-Pkm*H'*Kk'+Kk*Wk*Kk';
end




