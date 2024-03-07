%% housekeeping
clear all
close all
clc

%% get some calibration data
npts = 1000;
[t,yc] = serial_reader_n(npts,"COM5");%cm

%% do calibration stats
measuredDist = 15;%cm
v = var(yc);
meanDist = mean(yc)
offset = measuredDist - meanDist

%% get some test data
npts = 1000;
[t,y] = serial_reader_n(npts,"COM5");%cm


%% diff values
Q = [.01,.01,.05,.01,.01,.01]
R = [.1,.01,.1,.1,.1,v]
P0 = [0,0,0,10,0,0]
x0 = [0,0,0,0,10,0]

for i=1:length(R)
    runTest(y,t,Q(i),R(i),P0(i),x0(i),i)
end
%%
function runTest(y,t,Q,R,P0,x0,I)
% Kalman Inputs

F = 1;
G = 1;
H = 1;
Q = eye(1)*Q;
R = eye(1)*R;
L = [];

% P0 = 0
% x0 = y(1)

P = [P0];
xh = [x0];

u = 0;
u = u*ones(1,length(y));
% Kalman Filter

for k = 2:length(y)
    xhk_1 = xh(:,:,k-1);
    Pk_1 = P(:,:,k-1);
    [xhk,yhk,Pk] = Kfilter(xhk_1,Pk_1,u(k),y(k),F,G,H,Q,R);
    xh(:,:,k) = xhk;
    yh(:,:,k) = yhk;
    P(:,:,k) = Pk;
end

% plot kalman outputs
% plot(squeeze(yh))
figure
hold on 
plot(t,y)
plot(t,squeeze(xh))
legend(["Signal","Kalman Filter Estimation"])
xlabel("Sample Number (k)")
ylabel("Distance [cm]")
title("Q:"+string(Q)+" , "+"R:"+string(R)+" , "+"P0:"+string(P0)+" , "+"x0:"+string(x0))
end

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

function [t, data] = serial_reader_n(numDataPoints,portStr)

arguments
    numDataPoints
    portStr = 'COM9'
end

% Copyright 2014 The MathWorks, Inc.
% with some modificaitons by Eric Mehiel
% Create serial object for Arduino
if (~isempty(instrfind))
    fclose(instrfind); % close all ports to start, just to be sure...
    delete(instrfind);
end
s = serial(portStr); % change the COM Port number as needed
% Connect the serial port to Arduino

try
    fopen(s);
catch err
    fclose(instrfind);
    delete(instrfind);
    error('Make sure you select the correct COM Port where the Arduino is connected.');
end

% Read the data from Arduino
i = 0; % counter
%numDataPoints = 10;
data = zeros(numDataPoints,1);
t = zeros(numDataPoints,1);
tempData = fscanf(s); % Clear the serial buffer
timer = tic; % Start timer
while i < numDataPoints
    i = i + 1;
    % Read buffer data
    disp(string(i)+'/'+string(numDataPoints));
    data(i,1) = fscanf(s, '%f')'; % Change format string as needed
    % Read time stamp
    t(i) = toc(timer);
end
fclose(s);
delete(s);
clear s;
end
