clear all;
close all;
clc;

Cx = @(theta)...
    [[1 0 0];...
    [0 cosd(theta) sind(theta)];...
    [0 -sind(theta) cosd(theta)]];

Cy = @(theta)...
    [[cosd(theta) 0 -sind(theta)];...
    [ 0 1 0];...
    [sind(theta) 0 cosd(theta)]];

Cz = @(theta)...
    [[cosd(theta) sind(theta) 0];...
    [-sind(theta) cosd(theta) 0];...
    [0 0 1]];

C21 = Cx(180); % Ft -> Fa

% Fa0 will be ref frame a/c at takeoff      2
% Fa will be ref frame a/c at time t        3
% Ft will be ref frame ENU                  1

psi = 20; % z
theta = -5; % y
phi = 10; % x

% psi = 0; % z
% theta = 0; % y
% phi = 0; % x

%2a
C32 = Cx(phi)*Cy(theta)*Cz(psi); % Fa0 -> Fa

%2b
C21 = Cx(180);
C31 = C32*C21;

%2c
North1 = [0;1;0];
North3 = C31*North1;

%2d
eta =@(C) .5*sqrt(1+trace(C));
epsilon =@(C,eta)...
    [(C(2,3)-C(3,2))/(4*eta);...
    (C(3,1)-C(1,3))/(4*eta);...
    (C(1,2)-C(2,1))/(4*eta)];

C13=C31';

x=[1 0 0];
y=[0 1 0];
z=[0 0 1];

vecs(:,1)=x*C13;
vecs(:,2)=y*C13;
vecs(:,3)=z*C13;

eta = eta(C13);
vecs(:,4) = epsilon(C13,eta)';

Us = vecs(1,:);
Vs = vecs(2,:);
Ws = vecs(3,:);

figure
hold on
for i = 1:length(vecs)
    quiver3(0,0,0,Us(i),Vs(i),Ws(i));
end

legend("y","z","x","epsilon")
xlabel("E")
ylabel("N")
zlabel("U")
axis("equal")
title("Fa/c Orthonormal Basis and Epsilon In EMU")

%% Further Analysis
disp("The kinematic variables are defined within a given reference frame. To all other reference frames (including intermediate and transformed), they are relative. In this requard a kinematic variable is absolute in its reference frame until the reference frame begins to rotate.")



