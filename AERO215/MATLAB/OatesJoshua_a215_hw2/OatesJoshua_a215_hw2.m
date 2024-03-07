% Joshua Oates - a215 - hw2 - drag and power - fall 2021
clear all;
close all;
clc;

%constants
ftps_knots = 1.687810; %1 knot = 1.687810 ft/sec
gamma = 1.4; %(part 4)
R = 1716; %ft*lbf/(slug/rankin)
g = 32.2; %ft/s
%set up 4 cases

%baseline
baseline.WoS = 8.27;% Wing Loading W/S: 8.27 lbf/ft2
baseline.AR = 5.63;% Aspect Ratio AR: 5.63
baseline.CD0 = .0275;% Parasite Drag Coefficient CD,0: 0.0275
baseline.OE = .78;% Oswald Efficiency e: 0.78
baseline.CLmax = 1.4;% Clean Maximum Lift CL,max: 1.4
baseline.S = 127; % Wing Reference Area S: 127 ft2

%cruise
%change:
%nothing but run with different rho
cruise.WoS = 8.27;% Wing Loading W/S: 8.27 lbf/ft2
cruise.AR = 5.63;% Aspect Ratio AR: 5.63
cruise.CD0 = .0275;% Parasite Drag Coefficient CD,0: 0.0275
cruise.OE = .78;% Oswald Efficiency e: 0.78
cruise.CLmax = 1.4;% Clean Maximum Lift CL,max: 1.4
cruise.S = 127; % Wing Reference Area S: 127 ft2

%landing 
% Change:
% Parasite Drag Coefficient CD,0: 0.0975
% Maximum Lift CL,max: 1.7
landing.WoS = 8.27;% Wing Loading W/S: 8.27 lbf/ft2
landing.AR = 5.63;% Aspect Ratio AR: 5.63
landing.CD0 = .0975;% Parasite Drag Coefficient CD,0: 0.0975
landing.OE = .78;% Oswald Efficiency e: 0.78
landing.CLmax = 1.7;% Clean Maximum Lift CL,max: 1.7
landing.S = 127; % Wing Reference Area S: 127 ft2

%heavyload
% change
% Wing Loading (W/S): 10.39 lbf/ft2
heavyload.WoS = 10.39;% Wing Loading W/S: 10.39 lbf/ft2
heavyload.AR = 5.63;% Aspect Ratio AR: 5.63
heavyload.CD0 = .0275;% Parasite Drag Coefficient CD,0: 0.0275
heavyload.OE = .78;% Oswald Efficiency e: 0.78
heavyload.CLmax = 1.4;% Clean Maximum Lift CL,max: 1.4
heavyload.S = 127; % Wing Reference Area S: 127 ft2



%%Part2
%call function for hand calcs
HCspeed1knots = 70; %speed in knots for each senario
HCspeed1ftps = HCspeed1knots * ftps_knots;
sealevelrho = .00238;
cruiserho = .00189760; 
hold("on");
%for baseline
[Vhcb,Dhcb,Phcb]=dragPowerOatesJoshua(sealevelrho,HCspeed1ftps,baseline);
% plot(Vhcb,Dhcb);
%for cruise
[Vhcc,Dhcc,Phcc]=dragPowerOatesJoshua(cruiserho,HCspeed1ftps,cruise);
% plot(Vhcc,Dhcc);
%for landing
[Vhcl,Dhcl,Phcl]=dragPowerOatesJoshua(sealevelrho,HCspeed1ftps,landing);
% plot(Vhcl,Dhcl);
%for heavyload
[Vhch,Dhch,Phch]=dragPowerOatesJoshua(sealevelrho,HCspeed1ftps,heavyload);
% plot(Vhch,Dhch);
disp("Part 2:");
%disp(" ");
disp("The power at 70 knots for a baseline plane at sealevel is: " + Phcb(length(Phcb))+ " hp.");
disp("The drag  at 70 knots for a baseline plane at sealevel is: " + Dhcb(length(Dhcb))+ " lbs.");
%disp(" ");
disp("The power at 70 knots for a baseline plane at 7500ft is: " + Phcc(length(Phcc))+ " hp.");
disp("The drag at 70 knots for a baseline plane at 7500ft is: " + Dhcc(length(Dhcc))+ " lbs.");
%disp(" ");
disp("The power at 70 knots for a landing plane at sealevel is: " + Phcl(length(Phcl))+ " hp.");
disp("The drag at 70 knots for a landing plane at sealevel is: " + Dhcl(length(Dhcl))+ " lbs.");
%disp(" ");
disp("The power at 70 knots for a heavy plane at sealevel is: " + Phch(length(Phch))+ " hp.");
disp("The drag at 70 knots for a heavy plane at sealevel is: " + Dhch(length(Dhch))+ " lbs.");
disp(" ");

%%Part3
speed2knots = 140; %speed in knots for each senario
speed2ftps = speed2knots * ftps_knots;
sealevelrho = .00237717;
cruiserho = 0.00189760; %rho from digital dutch slug/ft^3
hold("on");
%for baseline
[V3b , D3b , P3b] =dragPowerOatesJoshua( sealevelrho, speed2ftps, baseline);
plot(V3b,D3b);
%for cruise
[V3c, D3c, P3c]=dragPowerOatesJoshua( cruiserho, speed2ftps, cruise);
plot(V3c , D3c);
%for landing
[V3l , D3l , P3l]=dragPowerOatesJoshua( sealevelrho, speed2ftps, landing);
plot(V3l,D3l);
%for heavyload
[V3h , D3h , P3h]=dragPowerOatesJoshua( sealevelrho, speed2ftps, heavyload);
plot(V3h,D3h);

V3b = V3b / ftps_knots;
V3c = V3c / ftps_knots;
V3l = V3l / ftps_knots;
V3h = V3h / ftps_knots;

figure(1)

subplot (1,2,1)
hold("on");
plot(V3b , D3b);
plot(V3c , D3c);
plot(V3l , D3l);
plot(V3h , D3h);
xlabel("Velocity (knots)");
ylabel("Drag (lbf)");
legend("base","cruise","landing","heavy");

subplot(1,2,2);
hold("on");
plot(V3b , P3b);
plot(V3c , P3c);
plot(V3l , P3l);
plot(V3h , P3h);
xlabel("Velocity (knots)");
ylabel("Power Required (HP)");
legend("base","cruise","landing","heavy");


%%Part4
%energy altitude
%In the same MATLAB script as used in Part 3, calculate the energy height of
% an aircraft  flying at 10,000 ft and Mach 0.7. (use gamma of 1.4). Calculate 
% the speed, in Mach number, the aircraft would be flying if it were to descend
% along a constant energy height contour from that initial speed and altitude 
% to sea-level. Calculate the starting and ending speed of sounds (at the two
% altitudes). Have MATLAB display all information asked for in question 4 in 
% the command window with units.

%calculate mach numbers
machi=.7; %mach number initial

Alti = 10000; %ft, altiltude initial
Altf = 0;%ft
Ti = 483;%Temp intial in rankin
Tf = 519;%rankin
ai = sqrt ( gamma * R * Ti ); %Speed of sound ft/s
af = sqrt( gamma * R * Tf );



%calculate energy height at alt 10000 and mach .7
Vi = machi * ai; %mach number times speed of sound, ft/s

%mechanical energy initial, sum of potential and kinetic:
zei = .5 * ( Vi ^ 2 ) / g + Alti ; % ft


%calculate energy height at alt 0 and mach (x)
zef = zei;
Vf = sqrt( ( zef - Altf ) * 2 * g );
machf = Vf / af;

%Output calculations
disp("Part 4:")
disp("The initial speed is mach 0.7 = "+Vi+" ft/s, at "+Alti+" ft, where the speed of sound is " + ai + " ft/s.")
disp("The final speed is mach " + machf + " = "+Vf+" ft/s, at "+Altf+" ft, where the speed of sound is " + af + " ft/s.")
disp("The energy height along this curve is " + zei + " ft.")



