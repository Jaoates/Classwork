%% final project mass properties solver
clear all
close all
clc


%% single cube

% detumble
detumbleMOI.m = 640; % mass in kg
detumbleMOI.dim = [2;2;2]; % x y and z lengths
detumbleMOI.CM = [0;0;0]; % x y and z locations of CM
detumbleMOI.name = "Detumble";
detumbleMOI.moi = moi(detumbleMOI);


%% individual parts

% bus
piece1.m = 500; % mass in kg
piece1.dim = [2;2;2]; % x y and z lengths
piece1.CM = [0;0;0]; % x y and z locations of CM
piece1.name = "Bus";

% sensor
piece2.m = 100; % mass in kg
piece2.dim = [.25;.25;1]; % x y and z lengths
piece2.CM = [0;0;1.5]; % x y and z locations of CM
piece2.name = "Sensor";

% panels
piece3.m = 20; % mass in kg
piece3.dim = [2;3;.05]; % x y and z lengths
piece3.CM = [0;2.5;0]; % x y and z locations of CM
piece3.name = "Panel1";

piece4 = piece3;
piece4.m = 20; % mass in kg
piece4.dim = [2;3;.05]; % x y and z lengths
piece4.CM = [0;-2.5;0]; % x y and z locations of CM
piece4.name = "Panel2";


pieces = [piece1,piece2,piece3,piece4];
clear piece1 piece2 piece3 piece4
%% do calcs

% moi pieces
for i = 1:length(pieces)
    pieces(i).moi = moi(pieces(i));
end

% CM nominal
nominal.m = sum([pieces.m]);
nominal.CM = sum([pieces.m].*[pieces.CM]./nominal.m,2);
nominal.name = "Nominal";


%% moi nominal - parallel axis

% distance from total CM to peice CM
for i = 1:length(pieces)
    r = nominal.CM - pieces(i).CM;

    x = r; % r is location of the peice CM in BOD
    x(1) = 0; % perpendicular location from x axis
    x = norm(x); % perpependicular distance from x axis

    y = r;
    y(2) = 0;
    y = norm(y);
    
    z = r;
    z(3) = 0;
    z = norm(z);

    pieces(i).r = [x;y;z]; % set of perpendicular distances from their respective axes
    
end

nominal.moi = sum([pieces.moi],2)+sum([pieces.m].*[pieces.r].^2,2);

%% fix Coords about CM of SC instead of geometric center of Bus
for i = 1:length(pieces)
    pieces(i).CM = pieces(i).CM - nominal.CM;
end

nominal.CM = nominal.CM-nominal.CM;

%% output

disp("The mass properties of the S/C in detumble operation are as follows: ")
disp("Mass in kg: "+string(detumbleMOI.m))
disp("The location of the ceter of mass in m in the body referenece frame in m: ")
disp(detumbleMOI.CM)
disp("Ix: "+string(detumbleMOI.moi(1)))
disp("Iy: "+string(detumbleMOI.moi(2)))
disp("Iz: "+string(detumbleMOI.moi(3)))

disp(" ")
disp(" ")

disp("The mass properties of the S/C in nominal operation are as follows: ")
disp("Mass in kg: "+string(nominal.m))
disp("The location of the ceter of mass in m in the body referenece frame in m: ")
disp(nominal.CM)
disp("The location of the geometric center of the bus in the body referenece frame in m: ")
disp(pieces(1).CM)
disp("Ix: "+string(nominal.moi(1)))
disp("Iy: "+string(nominal.moi(2)))
disp("Iz: "+string(nominal.moi(3)))

nominal.CenterOfBus = pieces(1).CM
nominal.magneticDipole = [0;0;1]*.5%A*m^2



% from nominal find other geo properties
a = 2; % cube side length
l_s = 1; % sensor length
a_s = .25; % sensor cube side length

l_p = 3; % solar panel length
w_p = 2; % width
t_p = .05; % thickness

A = [a^2 + 2*t_p*l_p + a_s*l_s, a^2 + a_s*l_s, a^2 + 2*l_p*w_p];
% projected area in body x, y, and z
nominal.A = [A A]; % same area in opposite directions

Ax = [a^2, 2*t_p*l_p, a_s*l_s];
z_b = [0, 0, 1.5];
rhoz = sum(Ax.*z_b)/sum(Ax);

nominal.rho = [[0; 0; rhoz]-nominal.CenterOfBus, [0; 0; rhoz]-nominal.CenterOfBus, [0; 0; 0]];
nominal.rho = [nominal.rho [[0; 0; rhoz]-nominal.CenterOfBus, [0; 0; rhoz]-nominal.CenterOfBus, [0; 0; 0]]];


%% Iwheel
Is = 1.2;
It = .6;
mw = 1;
nominal.Is = diag([Is,Is,Is])




save("MOInominal.mat","nominal","-mat")
save("MOIdetumble.mat","detumbleMOI","-mat")
%% functions
function out = moi(piece)
out(1) = (1/12)*piece.m*(piece.dim(2)^2+piece.dim(3)^2);
out(2) = (1/12)*piece.m*(piece.dim(1)^2+piece.dim(3)^2);
out(3) = (1/12)*piece.m*(piece.dim(1)^2+piece.dim(2)^2);
out = out';
end








