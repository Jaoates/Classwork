% Joshua Oates - a215 - Project 1 - Fall 2021
clear all;
close all;
clc;

%% common things, Planes structure setup
% Vans RV-12iS Parameters:
RV12.We = 740; % Empty Weight: 740 lbf
RV12.NumPass = 2; % Passengers: 2
RV12.Wf = 119; % Max Fuel load: 119 lbf
RV12.Altc = 7500; % Cruise altitude: 7,500 ft
RV12.CD0 = 0.0275;% Parasite Drag Coefficient CD0: 0.0275
RV12.AR = 5.63; % Aspect Ratio AR: 5.63
RV12.OE = 0.78; % Oswald Efficiency e: 0.78
RV12.CLmax = 1.4; % Clean Maximum Lift CLmax: 1.4
RV12.S = 127; % Wing Reference Area S: 127 ft2
RV12.Vmaxk =125; %125 knots

% Vans RV-10 Parameters:
RV10.We = 1605; % Empty Weight: 1605 lbf
RV10.NumPass = 4; % Passengers: 4
RV10.Wf = 360; % Max Fuel load: 360 lbf
RV10.Altc = 8000; % Cruise altitude: 8,000 ft
RV10.CD0 = 0.0265; % Parasite Drag Coefficient CD0: 0.0265
RV10.AR = 6.81; % Aspect Ratio AR: 6.81
RV10.OE = 0.79; % Oswald Efficiency e: 0.79
RV10.CLmax = 1.15;% Clean Maximum Lift CLmax: 1.15
RV10.S = 148; % Wing Reference Area S: 148 ft2
RV10.Vmaxk = 175; %175 knots

% Common Parameters:
PSFC_G_hrhp = .07;% Power specific fuel consumption 0.07 gallons per hour per horsepower
USDoGal = 5.04;% Avgas cost per gallon $5.04 (self serve at SLO airport)
WoGal = 6.01;% Avgas density 6.01 lbf per gallon
lbfoPass = 170;% 170 lbf weight per passenger

hpConv = 550; % 1 horsepower = 550 ft⋅lbf/sec
knotConv = 1.687810; % 1 knot = 1.687810 ft/sec
% 1 slug = 14.59390 kg
nMileConv = 6076; % 1 nautical miles = 6076 ft
PSFCConv = 6.01/(3600*550); % 1 gal/(hr*hp) = 6.01/(3600*550)1/ft

sealevelrho = 0.00238; %slug/ft^3 from dutch digital
c8000rho = 0.00186850; % RV-10
c7500rho = 0.00189760; % RV-12



%% Part 1

% estimate the minimum drag value (in lbf), the minimum drag speed (in knots),
% and the stall speed for the RV-12iS and the RV-10 at sea level at max takeoff weight
% (full fuel and full passengers) and then at sea level at min weight 
% (empty fuel,  full passengers). You can use MATLAB code to find the minimum
% values in the results or just take the minimum values off the plots.

% Create a figure with one subplot the sea level max weight drag comparison,
% one subplot of the sea level min weight drag comparison. Use 125 knots for the 
% RV-12iS max velocity and 175 knots for the RV-10 max velocity.

% Create a table with the values above and comment on the differences you notice.

%RV10 Sea Level full
% W is empty plus fuel plus two Pass 
RV10.WoS = (RV10.We + RV10.Wf +( lbfoPass * RV10.NumPass )) / RV10.S; 
[V10f,D10f,P10f,Mins10f,vStall10f] = dragPowerOatesJoshua (sealevelrho, RV10.Vmaxk * knotConv, RV10);

%RV12 Sea Level full
% W is empty plus fuel plus two Pass 
RV12.WoS = (RV12.We + RV12.Wf +( lbfoPass * RV12.NumPass )) / RV12.S; 
[V12f,D12f,P12f,Mins12f,vStall12f] = dragPowerOatesJoshua (sealevelrho, RV12.Vmaxk * knotConv, RV12);

%RV10 Sea Level empty
% W is empty plus two Pass 
RV10.WoS = (RV10.We +( lbfoPass * RV10.NumPass )) / RV10.S; 
[V10e,D10e,P10e,Mins10e,vStall10e] = dragPowerOatesJoshua (sealevelrho, RV10.Vmaxk * knotConv, RV10);

%RV12 Sea Level empty
% W is empty plus two Pass 
RV12.WoS = (RV12.We +( lbfoPass * RV12.NumPass )) / RV12.S; 
[V12e,D12e,P12e,Mins12e,vStall12e] = dragPowerOatesJoshua (sealevelrho, RV12.Vmaxk * knotConv, RV12);


subplot ( 1, 2, 1 );
hold("on");
plot ( V10f, D10f );
plot ( V12f, D12f );
xlabel("Velocity (knots)");
ylabel("Drag (lbf)");
legend("RV-10","RV-12");
title("at full / takeoff weight")

subplot (1 ,2 ,2 );
hold("on");
plot ( V10e, D10e );
plot ( V12e, D12e );
xlabel( "Velocity (knots)");
ylabel("Drag (lbf)");
legend("RV-10","RV-12");
title("at empty / landing weight")

%% Part 2

% The RV-12iS and the RV-10 have different nominal operational altitudes 
% ("cruise altitude") given above; estimate the minimum drag value (in lbf), 
% the minimum drag speed (in knots) for each aircraft at their respective 
% cruise altitudes, with full passenger load and half fuel.

% Add these values to the table from Part 1 and comment on the differences
% between the sea level and ceilings values.

% W is empty plus 1/2 fuel plus two Pass 
RV10.WoS = (RV10.We + (RV10.Wf / 2) +( lbfoPass * RV10.NumPass )) / RV10.S; 
[V10ch,D10ch,P10ch,Mins10ch,vStall10ch] = dragPowerOatesJoshua (c8000rho, RV10.Vmaxk * knotConv, RV10);

%RV12 Sea Level full
% W is empty plus 1/2 fuel plus two Pass 
RV12.WoS = (RV12.We + (RV12.Wf / 2) +( lbfoPass * RV12.NumPass )) / RV12.S; 
[V12ch,D12ch,P12ch,Mins12ch,vStall12ch] = dragPowerOatesJoshua (c7500rho, RV12.Vmaxk * knotConv, RV12);

%output to MinsXXch

%% Part 3
% Use the Breguet range equation to calculate range (in nautical miles) and 
% flight time (in hours) based on the half-full fuel and full passenger load.
% Remember, the Breguet range equation assumes constant velocity, lift-to-drag
% ratio, and specific fuel consumption (SFC). Therefore, assume steady level 
% flight and use the weight and minimum drag value from Part 2 to create an 
% estimate of lift-to-drag ratio during cruise. Use the minimum drag velocity
% from Part 2 as well. Use the max and min weights from Part 1 for the integration.

% With the range value, use the minimum drag velocity from Part 2 to estimate 
% the maximum flight time.

%TSFC = PSFC * V
%R = ( V / TSFC )*( L / D )* ln( Wi / Wf ) = PSFC *( L / D )* ln( Wi / Wf )
%L = Wh assmume SLF

PSFC = PSFC_G_hrhp * PSFCConv; %conversion to 1/ft
RV10.Wh = (RV10.We + (RV10.Wf / 2) +( lbfoPass * RV10.NumPass )); %weight at half fuel
RV12.Wh = (RV12.We + (RV12.Wf / 2) +( lbfoPass * RV12.NumPass ));
RV10.Wi = (RV10.We + (RV10.Wf) +( lbfoPass * RV10.NumPass )); %weight at takeoff
RV12.Wi = (RV12.We + (RV12.Wf) +( lbfoPass * RV12.NumPass ));
RV10.Wt = (RV10.We + ( lbfoPass * RV10.NumPass )); %weight at final
RV12.Wt = (RV12.We + ( lbfoPass * RV12.NumPass ));

RV10.R = (1/PSFC) * ( RV10.Wh / Mins10ch.D )* log( RV10.Wi / RV10.Wt ); %ft
RV12.R = (1/PSFC) * ( RV12.Wh / Mins12ch.D )* log( RV12.Wi / RV12.Wt ); 

RV10.R = RV10.R / nMileConv; 
RV12.R = RV12.R / nMileConv;

RV10.flightT = RV10.R / Mins10ch.V;
RV12.flightT = RV12.R / Mins12ch.V;
 
%% Part 4
% Using the cost of avgas given above and the results from Part 3, calculate
% the cost of the flight in total dollars spent on fuel,  then with that value
% find dollars per passenger, and dollars per passenger per nautical mile. 
% Assume all the fuel was consumed on the flight. 

% (Wf / WoGal) * USD/Gal = $
RV10.flightCost = (RV10.Wf / WoGal) * USDoGal;
RV12.flightCost = (RV12.Wf / WoGal) * USDoGal;

RV10.USDoPass = (RV10.flightCost / RV10.NumPass);
RV12.USDoPass = (RV12.flightCost / RV12.NumPass);

RV10.USDoPassNmile = (RV10.flightCost / (RV10.NumPass * RV10.R ));
RV12.USDoPassNmile = (RV12.flightCost / (RV12.NumPass * RV12.R ));






