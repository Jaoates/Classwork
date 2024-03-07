function [dv1,dv2,dv,T,ht,ecct,vt1,vt2] = joshHomann(r1,v1,r2,v2,mu)
% takes the magnitudes of r1 v1 at either apoapse or periapse of orbit 1
% and r2 v2 at the apoapse or periapse of orbit 2 as scalars. 
% it is assumed that the apses of orbits 1 and 2 are on oposite sides of 
% of the foci and that they lie on the same apse line.
% it is assumed that both orbits have the same grade, ie both pro- or retro-
% grade. parameters should be positive values but the returned dv's
% may be negative to corresponding to the retrograde burn. 
% The first 3 returned values correspond to delta V's
% the 4th returned value corresponds to transfer time ie 1/2 of period
% The 5th-8th returned values correspond to properties of the transfer orbit

arguments
    r1(1,1) double {mustBeNonnegative}
    v1(1,1) double {mustBeNonnegative}
    r2(1,1) double {mustBeNonnegative}
    v2(1,1) double {mustBeNonnegative}
    mu (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 3.986004418 * (10^5) %km^3/s^2
end
    ecct = ((r2-r1)/(r1+r2)); % absolute value so that if r2 > r1, ecct is postive
    ht = sqrt(mu*r1*(1+ecct)); % this assumes we're at periapse
    vt1 = ht/r1;
    vt2 = ht/r2;
    dv1 = vt1-v1;
    dv2 = v2-vt2;
    dv = abs(dv1)+abs(dv2);
    at = (r1+r2)/2;
    T = at^1.5*pi/sqrt(mu);
end

