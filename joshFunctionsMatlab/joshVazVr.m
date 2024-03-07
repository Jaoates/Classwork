function [vaz,vr,gamma] = joshVazVr(theta,ecc,h,mu)
% gives magnitude of azmuthal velocity and radial velocity
% takes theta ecc and h
% optionally takes mu for the center body
% assumes shperical body and 2 body
% angles in rad
arguments
    theta (1,1) double {mustBeReal,mustBeNonNan}
    ecc (1,1) double {mustBeReal, mustBeNonnegative}
    h (1,1) double {mustBeReal,mustBePositive}
    mu (1,1) double {mustBeReal,mustBeNonNan} = 3.986004418 * (10^5) %km^3/s^2 mu_earth
end
vr  = (mu/h)*ecc*sin(theta);
vaz = (mu/h)*(1+ecc*cos(theta));
gamma = atan2(vr,vaz);
end

