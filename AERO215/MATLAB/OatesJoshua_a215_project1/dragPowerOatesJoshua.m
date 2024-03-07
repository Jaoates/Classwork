function [V,D,P,minD,vStall] = dragPowerOatesJoshua(rho,vMax,inputs)
% calculate stall, drag, power for input aircraft
% velocity input in ft/sec (for some reason)
%rho in slug/ft^3
%returs in  knots, lbf, hp
%will return a strucutre with the minimum of D and corresponding V 
%as well for all three outputs.
% minD:
%   .V
%   .D


WoS=inputs.WoS;
AR=inputs.AR;
CD0=inputs.CD0;
OE=inputs.OE;
CLmax=inputs.CLmax;
S=inputs.S;

%calculate other values

W = WoS * S;

%rho = .00238; %slugs/ft^3


%calculate stall velocity
vStall = sqrt(W/(.5*rho*CLmax*S));

%from stall to Vmax calculate drag and power
V = linspace(vStall,vMax);
D = zeros(1,length(V));
P = zeros(1,length(V));

for i = 1:length(V)
    q = 0.5 * rho *V(i)^2;%dynamic pressure in lbf/ft^2
    Dp = CD0 * q * S; %parasitic drag in lbf
    Di = (W ^ 2) / ( q * S * pi * OE * AR ); %induced drag in lbf
    D(i) = Dp + Di; %total drag vector lbf
    P(i) = (D(i) * V(i))/550; %in hp, the 550 is the conversion factor        
end

V = V/1.687810; %convert fps to knots
[minD.D,I] = min(D);
minD.V = V(I);
    
    

