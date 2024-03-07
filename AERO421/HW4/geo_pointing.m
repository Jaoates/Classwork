function  [T,Te] = geo_pointing(r0, v0, T0lla, jd, err)
% function  [T,Te] = geo_monte_sim(r0, v0, p0hat, T0lla, jd, err)
% Te - T = error_vector
% T is actual target location
% Te is actual pointing location

% Mutate T0 -> T
T = T0lla;
T(1) = T(1)+err(2).value*randn; % lat
T(2) = T(2)+err(3).value*randn; % lon
T(3) = T(3)+err(4).value*randn; % alt
T = lla2eci_421(T(1),T(2),T(3),jd); % to ECI
T = T';

% recover T0 from T0lla
T0 = lla2eci_421(T0lla(1),T0lla(2),T0lla(3),jd + err(1).value*randn*1.15741e-5);
T0 = T0';

% get P0 from deffinition:
P0 = T0-r0;% note, P0 shouldn't really be used for anything, its just where the sc thinks its looking
P0hat = P0/norm(P0);

% Mutate r0 -> r
r = r0 + err(5).value*randn*joshRandSphere; % in track cross track radius (each has the same error)
r = r + err(8).value*randn*joshRandSphere; % sensor mounting error position (it actually just means r is off by a bit, but could also be applied to pointing)

% Mutate v0 -> v

% create temp basis for mutation with t3 in radial dircetion and t2 and t1 in
% the cross/in track plane
[t1,t2,t3] = joshFindPerpVec(r);

% get combination from 
temp = joshRandCircle;
v_error =randn*err(6).value*(temp(1)*t1 + temp(2)*t2); 
v = v0 + v_error;

% since t3 is in the radial direction, we just multiply the error by t3 and
% add that to v
v = v + err(7).value * randn * t3; % radial

% mutate r again using v and sidereal time error
% r = r + err(1).value*randn*v;
r = r + err(1).value*randn*v_error;


% this should be read as a rotation about a random axis that lies
% perpendicular to P0 so ill get that axis the same was as from before
[t1,t2] = joshFindPerpVec(r0);
temp = joshRandCircle;
a = temp(1)*t1 + temp(2)*t2;
phi = deg2rad(err(9).value)*randn; % angle of rotation in rad
C = joshPrincAxe2RotM(a,phi);
Phat = C*P0hat;

Re = norm(T);
P = joshSphereIntersect(Phat,-r,Re);

Te = r+P;

end