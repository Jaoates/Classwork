function [] = dispCOEsOatesJoshua(COEs ,periodBool ,epsilonBool )
%Takes COEs Object and creates formated output. 
% Must have these properties in units of km , km/s , s , degrees
% COEs.a = 0;
% COEs.e = 0;
% COEs.nu = 0;
% COEs.i = 0;
% COEs.raan = 0;
% COEs.aop = 0;
% COEs.T = 0;
arguments
    COEs;
    periodBool logical = 0;
    epsilonBool logical = 0;
end

disp("------- " + COEs.name + " -------")
disp("Semi-major axis       : " + COEs.a + " km")
disp("eccentricity          : " + COEs.e)
disp("true anomaly          : " + COEs.nu + " degrees")
disp("inclination           : " + COEs.i + " degrees")
disp("RAAN                  : " + COEs.raan + " degrees")
disp("Argument of periapsis : " + COEs.aop + " degrees")

if periodBool == true
    disp("Period                : " + COEs.T + " seconds")
end

if epsilonBool == true
    disp("Specific Mech. Energy : " + COEs.E + " km^2/s^2")
end

disp(" ")

end

