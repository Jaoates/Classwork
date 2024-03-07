function [TRatio,PRatio,rhoRatio] = joshIsentropicToolbox(M,gamma)
% takes Mach number where M = velocity/speed of sound
% takes gamma as a ratio of Cp/Cv
% returns ratios of T0 P0 and rho0 over T P and rho where the numerator is T0 and the denominator is T 

arguments
    M (1,1) {mustBeNonnegative}
    gamma (1,1) {mustBeNonnegative} =1.4 % default gamma for air
end
TRatio = 1+((gamma-1)/2)*M^2;
PRatio = TRatio^(gamma/(gamma-1));
rhoRatio = PRatio*TRatio^-1;
end