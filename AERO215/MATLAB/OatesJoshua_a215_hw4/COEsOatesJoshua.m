function [a,em,nu,i,raan,aop,T,E] = COEsOatesJoshua(R,V,u)
%COESOATESJOSHUA Takes postion and velocity vector and returns COES, all in
%ECI fram of reffrenece, km and seconds as units and degrees
% a = semi major axis, e = eccentricity , i = inclination
%  raan = right accention acending node, aop = argument of periapsis
% will return T in s as a period if given 7 spots for outputs 
% and E in km^2/s^2 as specific mechanical energy if given 8



    %argName1 (dimensions) dataType {validators} = defaultValue
    % x (1,:) {mustBeNumeric}
    % defualt u is u for earth

    arguments
        R (1,3) {mustBeNumeric};
        V (1,3) {mustBeNumeric};
        u (1,1) {mustBeNumeric} = 3.986004418 * (10^5); %km^3/s^2
    end

% uearth = 3.986004418(8)	× 10^14 m^3/s^2
ihat=[1,0,0];
jhat=[0,1,0];
khat=[0,0,1];

Rm = norm(R);
Vm = norm(V);

%calculate orbital constants
h = cross(R,V); %angular momentum vector
hm = norm(h);
E = ( ( Vm^2 ) / 2) - ( u / Rm ) ; %specfic mechanical energy
%calculate COEs

a = -u / (2 * E ); %semi major axis in km
T = 2*pi*sqrt((a^3)/u); %period in s

e = (1/u) * (((Vm^2)-(u/Rm) ) * R - (dot(R,V) * V)); %eccentricity vector
em = norm (e); %magnitude of e

i = acosd((dot(khat,h))/hm); %inclination degrees
n = cross(khat,h); %node vector
nm = norm(n); %magnitude n

%raan
raan = acosd(dot(ihat,n)/nm);
if n(2) < 0  %checks the vector relative to j to see if angle is positive or negative
    raan = 360 - raan;
end

%aop
aop = acosd(dot(n,e)/(nm*em));
if e(3) < 0
    aop = 360 - aop;
end

%nu
nu = acosd(dot(e,R)/(em*Rm));
if(dot(R,V) < 0); % cehck flight path angle to see if it is postive or negative
    nu = 360 -nu;
end




end

