% Aero302: Aerospace Fluid Mechanics
% Lab 4: Water Tunnel
% Fall 2022

clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Colton Crosby and Justin Self %%%%
%%% Cal Poly Aerospace Engineering %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Frequency Data from Experiment

% Colton Crosby // Justin Self 
% 11/15/2022

% constants and equations
FPS = 30; % frames per second
rho = 1000; % kg/m3 water
mu = 0.001; % mu water Pa.s
MotorSpeed = 20; % Hz
u = 1.85/39.37; % inch/s to m/s (from motor speed calibration curve)

St = @(f,D) (f*D)/u; % strouhal number
Re = @(D) (rho*u*D)/mu; % reynolds number

% diameters
cyl4.D  = (1/1000)*7.82;  % m
cyl5.D  = (1/1000)*15.88; % m
cyl6.D  = (1/1000)*21.16; % m
cyl7.D  = (1/1000)*26.67; % m
cyl8.D  = (1/1000)*48.61; % m
cyl9.D  = (1/1000)*88.88; % m
cyl10.D = (1/1000)*114.42; % m

%% Cylinder 4  ----------------------------------------

% freqs & standard deviation   
cyl4.G.freq = 1/3.4625;                 % Hz
cyl4.R.freq = 1/3.527083333;            % Hz

cyl4.G.stdDev = 0.304170744;                      
cyl4.R.stdDev = 0.270039435;

%--------- Re and St
% green
cyl4.G.St = St(cyl4.G.freq, cyl4.D);    % Strouhal #
cyl4.G.Re = Re(cyl4.D);                 % Reynolds #

% red
cyl4.R.St = St(cyl4.R.freq, cyl4.D);    % Strouhal #
cyl4.R.Re = Re(cyl4.D);                 % Reynolds #

%% Cylinder 5  ----------------------------------------

% freqs & standard deviation               
cyl5.G.freq = 1/5.641666667;            % Hz 
cyl5.R.freq = 1/5.416666667;            % Hz

cyl5.G.stdDev = 0.489822341;
cyl5.R.stdDev = 0.290338245;

% Green dye
cyl5.G.St = St(cyl5.G.freq, cyl5.D);    % Strouhal #
cyl5.G.Re = Re(cyl5.D);                 % Reynolds #

% Red dye
cyl5.R.St = St(cyl5.R.freq, cyl5.D);    % Strouhal #
cyl5.R.Re = Re(cyl5.D);                 % Reynolds #

%% Cylinder 6  ---------------------------------------- 

% freqs & standard deviation               
cyl6.G.freq = 1/8.04;                   % Hz 
cyl6.R.freq = 1/6.696296296;            % Hz        

cyl6.G.stdDev = 0.641785703;
cyl6.R.stdDev = 0.54631199;

% Green dye
cyl6.G.St = St(cyl6.G.freq, cyl6.D);    % Strouhal #
cyl6.G.Re = Re(cyl6.D);                 % Reynolds #

% Red dye
cyl6.R.St = St(cyl6.R.freq, cyl6.D);    % Strouhal #
cyl6.R.Re = Re(cyl6.D);                 % Reynolds #

%% Cylinder 7  ---------------------------------------- 

% freqs & standard deviation               
cyl7.G.freq = 1/8.586666667;            % Hz 
cyl7.R.freq = 1/7.491666667;            % Hz        

cyl7.G.stdDev = 0.373124942;
cyl7.R.stdDev = 0.39569442;

% Green dye
cyl7.G.St = St(cyl7.G.freq, cyl7.D);    % Strouhal #
cyl7.G.Re = Re(cyl7.D);                 % Reynolds #

% Red dye
cyl7.R.St = St(cyl7.R.freq, cyl7.D);    % Strouhal #
cyl7.R.Re = Re(cyl7.D);                 % Reynolds #

%% Cylinder 8 (Colton Video 1098, first cyl (WHITE))  ---------------------------------------- 

% freqs & standard deviation               
cyl8.G.freq = 1/13.2;                   % Hz 
cyl8.R.freq = 1/10.08333333;            % Hz        

cyl8.G.stdDev = 0.670544279;
cyl8.R.stdDev = 0.535723809;

% Green dye
cyl8.G.St = St(cyl8.G.freq, cyl8.D);    % Strouhal #
cyl8.G.Re = Re(cyl8.D);                 % Reynolds #

% Red dye
cyl8.R.St = St(cyl8.R.freq, cyl8.D);    % Strouhal #
cyl8.R.Re = Re(cyl8.D);                 % Reynolds #

%% Cylinder 9 (Colton Video 1098, second cyl (BLACK))  ---------------------------------------- 

% freqs & standard deviation               
cyl9.G.freq = 1/10.61111111;            % Hz 
cyl9.R.freq = 1/11.32222222;            % Hz        

cyl9.G.stdDev = 0.435039377;
cyl9.R.stdDev = 1.124145766;     % (!) only 3 data points (!)

% Green dye
cyl9.G.St = St(cyl9.G.freq, cyl9.D);    % Strouhal #
cyl9.G.Re = Re(cyl9.D);                 % Reynolds #

% Red dye
cyl9.R.St = St(cyl9.R.freq, cyl9.D);    % Strouhal #
cyl9.R.Re = Re(cyl9.D);                 % Reynolds #

% % % % freqs & standard deviation              
% % % cyl9.G.freq = 1/9.533333333;            % Hzcyl9.R.freq = 1/11.32222222;            % Hz       
% % % cyl9.G.stdDev = 0.427741701;
% % % cyl9.R.stdDev = 1.124145766;     % (!) only 3 data points (!)
% % % 
% % % % Green dye
% % % cyl9.G.St = St(cyl9.G.freq, cyl9.D);    % Strouhal #
% % % cyl9.G.Re = Re(cyl9.D);                 % Reynolds #
% % % 
% % % % Red dye
% % % % cyl9.R.St = St(cyl9.R.freq, cyl9.D);    % Strouhal #
% % % % cyl9.R.Re = Re(cyl9.D);                 % Reynolds #
% % %  

%% Cylinder 10 (Colton Video 1098, third cyl (white; massive))  ---------------------------------------- 

% freqs & standard deviation               
cyl10.G.freq = 1/15.27777778;            % Hz 
cyl10.R.freq = 1/14.78888889;            % Hz        

cyl10.G.stdDev = 0.323751162;
cyl10.R.stdDev = 0.269430126;

% Green dye
cyl10.G.St = St(cyl10.G.freq, cyl10.D);    % Strouhal #
cyl10.G.Re = Re(cyl10.D);                 % Reynolds #

% Red dye
cyl10.R.St = St(cyl10.R.freq, cyl10.D);    % Strouhal #
cyl10.R.Re = Re(cyl10.D);                 % Reynolds #

%% Re, St plots

DM = [cyl4,cyl5,cyl6,cyl7,cyl8,cyl9,cyl10];
n=length(DM);
Rel = zeros(1,n);
Stl = zeros(1,n);

for i=1:n
    Rel(i) = DM(i).G.Re;
    Stl(i) = DM(i).G.St;
    disp(DM(i).G.stdDev)
end
Re = [Rel(1:5),Rel(7)];% remove outlyer
St = [Stl(1:5),Stl(7)];


p = polyfit(Re,St,1);
X = linspace(Re(1)-100,Re(end)+100,100);
Y = polyval(p,X);
SS = @(d) sum(d.^2);
res = polyval(p,Re)-St;
tot = polyval(p,Re)-mean(St);
R2 = 1-SS(res)/SS(tot); % R^2 of best fit line


figure
hold on
plot(X,Y)

for i=1:n-1
  scatter(Re(i),St(i),'filled')  
end
scatter(Rel(6),St(6),'x','red') % plot outlyer

xlabel("Reynolds Number",'Interpreter','latex')
ylabel("Strouhal Number",'Interpreter','latex')
str = [4,5,6,7,8,10];
str = string(str);
str = ["Best Fit",str,"9 (Outlier)"];
legend(str,'Interpreter','latex','Location','northwest')



