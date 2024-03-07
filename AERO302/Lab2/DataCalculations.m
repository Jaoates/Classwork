5%% 0 
clear all;
close all;
clc

%% Imports
SV1 = importdata("cylindertest1.mat"); % Scanivalve
SV2 = importdata("cylindertest2.mat"); 
SV3 = importdata("cylindertest3.mat");

SV4 = importdata("cylindertest4.mat");
SV5 = importdata("cylindertest5.mat");
SV6 = importdata("cylindertest6.mat");

V_sc = importdata("V_sc.mat");

DM = [SV1 SV2 SV3 SV4 SV5 SV6];

clear SV1 SV2 SV3 SV4 SV5 SV6
rho = 1.204; % kg/m^3

% from lab view the first port is tot second is stat, 3-26 is cylinder

rho = 1.204; % kg/m^3

for i = 1:length(DM)
    DM(i).V_sc = V_sc(i); % speed controller
    DM(i).V_sv = DM(i).V; % scanivalve
    DM(i).V_sv = mean(DM(i).V_sv);
    DM(i).P_gauge = mean(DM(i).P); % ?
    DM(i).P = DM(i).P_gauge + DM(i).Amb_Press;
    DM(i).P_tot = DM(i).P(1); % total pressure at test section
    DM(i).P_sta = DM(i).P(2); % static pressure at test section
    DM(i).P_dyn = DM(i).P_tot-DM(i).P_sta; % dynamic pressure
    DM(i).P_gauge = DM(i).P_gauge(3:end);
    DM(i).P = DM(i).P(3:end);
end

RPM = [500,600,300,500,600,300];
RPM = string(RPM);
for i = 1:length(DM)
    DM(i).RPM = RPM(i);
end
clear RPM
DM = rmfield(DM,'V');

for i = 1:length(DM)
    DM(i).V = sqrt((DM(i).P_dyn*2)/rho); % V at inlet from dynamic pressure
end

DM = rmfield(DM,"filepath");
DM = rmfield(DM,"t");

clear V_sc

% 
% figure
% Y = DM(5).P_gauge;
% bar(Y)
% figure
% Y = DM(2).P_gauge;
% bar(Y)

%% setting up port angles
n = 24;
ang_inc = (2*pi)/24;
prt_ang=zeros(1,n); % port angle
prt_unt=zeros(2,n); % port unit vector
for i = 1:n
    prt_ang(i) = ang_inc*(i-1); % here index 1 goes with port 24 since port 24 is nominally on the x axis
end 
prt_ang = prt_ang+ang_inc; % now the angles can map directly to the 
prt_ang(n) = 0;

for i = 1:n
    prt_unt(:,i) = -[cos(prt_ang(i));sin(prt_ang(i))]; % inward facing vectors that will be in the direction of force
end

dia = 155; % diameter in mm
dia = dia/1000;
cir = pi*dia; % circumfrenc in m
len = cir/24; % length of each section
A = 1*len; % area of each panel m^2, asumme cylinder is 1m long
clear  cir len
for i = 1:length(DM)
    DM(i).prt_frc=zeros(2,n); % force vector
    P1 = DM(i).P;
    P2 = P1;
    P2(23) = (P2(22)+P2(24))/2; % 23 is an outlyer
    DM(i).P2 = P2;
    
    for j =1:n
        DM(i).prt_frc1(:,j) = prt_unt(:,j).*P1(j)*A; % force vector acting per unit area on each port
        DM(i).prt_frc2(:,j) = prt_unt(:,j).*P2(j)*A;
    end
    DM(i).res_frc1 = [sum(DM(i).prt_frc1(1,:));sum(DM(i).prt_frc1(2,:))];
    DM(i).res_frc2 = [sum(DM(i).prt_frc2(1,:));sum(DM(i).prt_frc2(2,:))];
end

clear A P1 P2 i j ang_inc %prt_unt prt_ang

%% Cp vs theta
for i = 1:length(DM)
    P_inf = DM(i).P_sta; % freestream absolute static pressure, an averaged value from station 2 
    V = DM(i).V; % freestream Velocity
    CP = @(P) (P-P_inf)/(.5*rho*V^2); % definition of CP(P) where P is absolute pressure measured at a given port
    for j = 1:n
        DM(i).CP(j) = CP(DM(i).P(j)); % calculate CP for all ports
    end        
end

prt_ang(n) = 2*pi;

figure
hold on
for i = 1:length(DM)
    plot(prt_ang,DM(i).CP)
end
legend(DM.RPM)

%%%%%%%%%%%redo with P2

for i = 1:length(DM)
    P_inf = DM(i).P_sta;
    V = DM(i).V;
    CP = @(P) (P-P_inf)/(.5*rho*V^2);
    for j = 1:n
        DM(i).CP(j) = CP(DM(i).P2(j));
    end        
end

prt_ang(n) = 2*pi;

figure
hold on
for i = 1:length(DM)
    plot(prt_ang,DM(i).CP)
end
legend(DM.RPM)

%% Re vs cd
T = 20; % C

muRef = 1.716e-5;
TRef = 273.15;
S = 110.4;
mu = @(T) muRef*(((T+S).^-1).*(TRef+S)).*(T/TRef).^1.5;
T=T+TRef;
mu = mu(T);

A = dia;
Reynolds = @(dat) (rho*dia*dat.V)/mu;
Cd = @(dat) dat.res_frc2(1)/(.5*rho*A*dat.V^2);

figure
hold on
for i = 1:length(DM)
    DM(i).Reynolds = Reynolds(DM(i));
    DM(i).Cd = Cd(DM(i));
    X(i) = DM(i).Reynolds;
    Y(i) =  DM(i).Cd;
end



plot([X(3),X(1),X(2)],[Y(3),Y(1),Y(2)],"-o")
plot([X(6),X(4),X(5)],[Y(6),Y(4),Y(5)],"-*")
xlabel('Re','Interpreter','latex')
ylabel('$C_{D}$','Interpreter','latex')
legend("Sym","Trip")

%% quiver of resultants

vecs(:,1)=DM(2).res_frc2
vecs(:,2)=DM(5).res_frc2
vecs(:,3)=[10,0]

Us = vecs(1,:);
Vs = vecs(2,:);

figure
hold on
for i = 1:length(vecs)
    quiver(0,0,Us(i),Vs(i));
end

legend("Sym = <"+DM(2).res_frc2(1)+","+DM(2).res_frc2(2)+"> N","Trip = <"+DM(5).res_frc2(1)+","+DM(5).res_frc2(2)+"> N","Freestream Direction",'Location','best')
xlabel("Drag Force [N]")
ylabel("Lift Force [N]")
% myString = ("Sym = <"+DM(2).res_frc2(1)+","+DM(2).res_frc2(2)+"> N"+newline+"Trip = <"+DM(5).res_frc2(1)+","+DM(5).res_frc2(2)+"> N")
% dim = [.2 .7 .35 .1]
% annotation('textbox',dim,'String',myString)

%% Functions
function [av] = myav(vec)
    vec = sort(vec);
    vec(end) = [];
    vec(1) = [];
    av = mean(vec);
end