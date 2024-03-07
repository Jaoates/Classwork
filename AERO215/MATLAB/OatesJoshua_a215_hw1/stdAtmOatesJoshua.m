function [T,P,rho] = stdAtmOatesJoshua(h)
% returns various properties of the atmosphere for a given altitude
%   calculate Temperature in K, Pressure in Pa, and Desnsity in kg/m^3, for an altitude in m,
%   uses the 1959 model of the atmosphere


% Sea level temperature 288.16 K
% Sea level pressure 101.325 kPa
% Sea level density 1.2250 kg/m^3
%  
% 
% Gradient layer  -0.0065 K/m up to 11000m
% Isothermal layer  11000m to 25000m
% Gradient layer  0.0030 K/m up to 47000m
% Isothermal layer  47000m to 53000m
% Gradient layer  -0.0045 K/m up to 79000m
% Isothermal layer  79000m to 90000m
% Gradient layer  0.0040 K/m up to 100000m


%set parameters at sea level

T0 = 288.16; %K
P0 = 101325; %Pa
rho0 = 1.2250; %kg/m^3
h0 = 0; %m
g0 = 9.80665; %m/s2
R = 287; %J/(kg*K) 


%lapse rate per layer and top of layer
% a1 =  -0.0065; %K/m
% h1 = 11000; %m
% 
% a2 = 0; %K/m
% h2 = 25000; %m
% 
% a3 = 0.0030; %K/m
% h3 = 47000; %m
% 
% a4 = 0; %K/m
% h4 = 53000;%m
% 
% a5 = -0.0045; %K/m
% h5 = 79000;%m
% 
% a6 = 0; %K/m
% h6 = 90000;%m
% 
% a7 = 0.0040;  %K/m
% h7 = 100000;%m
% %solve for T, P, rho at each base

aMat = [-0.0065,0,0.0030,0,-0.0045,0,0.0040]; %a 1-7 in K/m 
hMat = [0,11000,25000,47000,53000,79000,90000,100000]; %h 0-7 in m

TMat = [T0];       %temp at base of layer 1-7 and at the top of 7
PMat = [P0];       %pressure at base of layer 1-7 and at the top of 7
rhoMat = [rho0];   %density at base of layer 1-7 and at the top of 7


for i = 1:7 % for each layer 1-7 calculate T P and rho
    TMat(i+1) = TMat(i) + aMat(i) * (hMat(i+1) - hMat(i));
    
    if aMat(i)==0 %isothermal
        PMat(i+1) = PMat(i)*(exp(-(g0/(R*TMat(i+1)))*(hMat(i+1)-hMat(i))));
        rhoMat(i+1) = rhoMat(i)*(exp(-(g0/(R*TMat(i+1)))*(hMat(i+1)-hMat(i))));
    else%gradient
        PMat(i+1) = PMat(i)*((TMat(i+1)/TMat(i))^(-g0/(aMat(i)*R)));
        rhoMat(i+1) = rhoMat(i)*((TMat(i+1)/TMat(i))^((-g0/(aMat(i)*R)-1)));
    end
end


inLayer = 0; %inLayer will hold the layer that that h is in


while hMat(inLayer+1) < h
    inLayer=inLayer+1;
end

if h == 0
    layer = 1;
end



%solve T, P, rho at h by resusing the formulas above. Now the output is
%final, use (inLayer) where i was before, use h T once solved

T = TMat(inLayer) + aMat(inLayer) * (h - hMat(inLayer));

if aMat(inLayer)==0 %isothermal
    P = PMat(inLayer)*(exp(-(g0/(R*T))*(h - hMat(inLayer))));
    rho = rhoMat(inLayer)*(exp(-(g0/(R*T))*(h - hMat(inLayer))));
else %gradient
    P = PMat(inLayer)*((T/TMat(inLayer))^(-g0/(aMat(inLayer)*R)));
    rho = rhoMat(inLayer)*((T/TMat(inLayer))^((-g0/(aMat(inLayer)*R)-1)));
end

end

