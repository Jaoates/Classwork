function [r,v] = joshCOE2rv(a,ecc,theta,inc,raan,aop,mu)
% all angles in rads

warning("Something seems to be wrong here run check cases before using")

[Cx,Cy,Cz]=joshAxisRotation(); % just rotation matrix about x and z of theta

rp = a*(1-ecc^2)/(1+ecc); % cos(0) = 1
h = sqrt(mu*rp*(1+ecc*cos(theta)));

r = (h^2/mu)/(1+ecc*cos(theta));
rperi = [cos(theta);sin(theta);0]*r; % r vector in perifocal
[vaz,vr] = joshVazVr(theta,ecc,h,mu); % v vector in local horizontal vertical
vloc = [vr;vaz;0];

vperi = Cz(-theta)*vloc; % v vector rotated to perifocal

Q = Cz(raan)*Cx(inc)*Cz(aop); 
Q = Q'; % Perifocal -> ECI

r = Q*rperi; % rotation to get to Inertial frame
v = Q*vperi;

r = r';
v = v';
end