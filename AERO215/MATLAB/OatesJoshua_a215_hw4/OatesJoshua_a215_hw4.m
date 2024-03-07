%% Joshua Oates - a215 - Fall 2021 - HW 4  - changing orbits
clear all;
close all;
clc;
u = 3.986004418 * (10^5); %km^3/s^2, (mu) earth

%% Part 1
%set initial R and V vectors
%call COEs function and display output
% Initial Position (radius): [6161.56, 454.07, -2537.72] km
% Initial Velocity: [0.376, 7.391, 2.224] km/s
% these vars will be erased when they are moved into a COEs struct

Rv = [6161.56, 454.07, -2537.72];
Vv = [0.376, 7.391, 2.224];


%create a template for a structure that can hold the COEs for a given orbit
COEs.name = "name0";
COEs.a = 0;
COEs.e = 0;
COEs.nu = 0;
COEs.i = 0;
COEs.raan = 0;
COEs.aop = 0;
COEs.T = 0;
COEs.E = 0;
COEs.Rv = [0,0,0];
COEs.Vv = [0,0,0];


COEs1 = COEs;% use template to create COEs1. COEs1 is for the initial orbit
COEs1.name = "COEs1: INITIAL ORBIT";
COEs1.Rv = Rv;
COEs1.Vv = Vv;
clear Rv;
clear Vv;

[COEs1.a,COEs1.e,COEs1.nu,COEs1.i,COEs1.raan,COEs1.aop,COEs1.T,COEs1.E] = COEsOatesJoshua(COEs1.Rv,COEs1.Vv);
dispCOEsOatesJoshua(COEs1,1,1);


%% Part 2
% make orbit geostat
% geostationary orbit (position magnitude of 42157 km and velocity magnitude of 3.07 km/s).
% solve for semi major axis
% solve for specific mechanical energy
% assume there is no k-hat component for final V or R and make one have
% only i-hat and the other j-hat ex. V=[0,v,0] R=[r,0,0]
COEs2 = COEs;
COEs2.name = "COEs2: FINAL ORBIT";
COEs2.Rv = [42157 ,0 ,0];
COEs2.Vv = [0, 3.07 ,0];

% calculate epsilon = specific mechanical energy = E
[COEs2.a,COEs2.e,COEs2.nu,COEs2.i,COEs2.raan,COEs2.aop,COEs2.T,COEs2.E] = COEsOatesJoshua(COEs2.Rv,COEs2.Vv);

dispCOEsOatesJoshua(COEs2,1,1);

%% Part 3
% Use 4 - burn method (circularize at apogee, simple plane change: i=0, hohmann transfer)
% Calculate the dV required for this method of transfer
% use epsilon equation to move between values.

% First Burn - circularize
% Find Ra for Orbit 1
COEs1.Ra = COEs1.a * ( 1 + COEs1.e); %km, Radius as magnitude at apoapsis

% E already calculated

% Find VRam - Velocity magnitude at apoapsis for Orbit 1
COEs1.VRa = sqrt( ( COEs1.E + ( u/COEs1.Ra) ) *2 );

% Find E for an orbit where Ra = COEs1.Ra and e = 0. ie circularize

COEsC.Ra = COEs1.Ra; %Orbit Circular ( will be left with actual COEs uncalculated... its a bit of a misnomer I guess )

% Find E of circular
COEsC.E = (-u) / ( 2* COEsC.Ra );

% Use E to find Vm ( magnitude V ) for circular orbit ( Orbit C )
COEsC.VRa = sqrt( ( COEsC.E + ( u/COEsC.Ra) ) *2 );

% Use V of orbit 1 and V of orbit C to find dV of burn 1 for method 1
dV1.b1 = abs( COEsC.VRa - norm(COEs1.Vv)  ); %km/s

% Burn 2 - simple plane change
% Set inclination to 0
COEsE.i = 0;
thetaCE = ( COEs1.i - COEsE.i ); %degrees (delta i?)
% calculate dv for burn 2 
dV1.b2 = 2 * COEsC.VRa * sind( thetaCE/2 );

% Burn 3 - HT1
% get a of transfer orbit ( orbit T )
COEsT.a = ( COEsC.Ra + COEs2.a )/2; 

% Find E for orbit T
COEsT.E = (-u) / ( 2* COEsT.a );

% Find velocity of T at Rp ( intersection with Orbit E ) as magnitude = VRp
COEsT.VRp = sqrt( ( COEsT.E + ( u/COEsC.Ra) ) *2 );

% Use VRa of Circular and VRp of transfer to find dV burn 3
dV1.b3 = abs( COEsT.VRp - COEsC.VRa );

% Burn 4 - HT2
% Use E to find VRa of orbit T
COEsT.VRa = sqrt( ( COEsT.E + ( u/norm(COEs2.Rv)) ) *2 );

%Find dV 
dV1.b4 = abs( norm(COEs2.Vv)- COEsT.VRa );

dV1.net = (dV1.b1 + dV1.b2 + dV1.b3 + dV1.b4 );

disp("In units of km/s , dV's for - circularizatoin - plane change - HT1 - HT2 - are :");
disp(dV1);
disp(" ");

%% Part 4
% combined burn method of transfer - circ, combined plane change and HT1, HT2 

% I know that nothing is different until the combined PC_HT1 so I will
% start at Orbit C again. I will also be going to the same Orbit T while
% skipping Orbit E. Orbit 2 = geo is also the same so, I only need to
% calculate PC_HT1 and I will add that dV to the previos dVs for HT2 and
% Circ.



% Using Law of cosines, we know that dV for PC_HT1 is 
% dV2.b2 = sqrt( Vi^2 +Vf^2 - ( 2* Vi * Vf * cos(Theta) ) )

dV2.b1 = dV1.b1; %circularize is same

dV2.b2 = sqrt( (COEsC.VRa ^ 2 + COEsT.VRp ^ 2) - ( 2* COEsC.VRa * COEsT.VRp * cosd(thetaCE) ) );

dV2.b3 = dV1.b4; %HT2 is same
dV2.net = ( dV2.b1 + dV2.b2 + dV2.b3 );

disp("In units of km/s , dV's for - circularizatoin - PC / HT1 - HT2 - are :");
disp(dV2);
