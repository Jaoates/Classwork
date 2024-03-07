
function [fz,y,A,z,flag,glag,gdotlag] = joshfLambert(r1,r2,dt,theta,mu,C,S)%,z0)%,pr)
%{
this function will return the lagrange coeffs along with z from the
universal variable, A and yfrom lamberts problem and fz from lamberts
problem such that the zero of fz will give you the solution z to lamberts
problem. r1 r2 should be scalars, dt is transit time
%}


arguments
    r1 (1,1) double {mustBeReal, mustBePositive}
    r2 (1,1) double {mustBeReal, mustBePositive}
    dt (1,1) double {mustBePositive}
    theta (1,1) double {mustBePositive}
    mu (1,1) double {mustBePositive} = 398600
    C (1,:) double {mustBeReal} = nan
    S (1,:) double {mustBeReal} = nan
    %     z0 (1,1) double {mustBeReal} = nan
end


% warning("joshfLambert: This function may be useful but it is not well tested and complete argument validation has not been implimented.")


if isnan(C)|isnan(S)
    [C,S] = joshStumpffCoeffs();
end

coefs = length(C);
if length(S)~=coefs
    throw(MException("joshfLambert:invalidInput","S and C should be the same length"))
end

A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));

C = @(z) sum(C.*joshStumpffZ(z,coefs));
S = @(z) sum(S.*joshStumpffZ(z,coefs));

y = @(z)(r1+r2+  A*((z*S(z)-1))/sqrt(C(z))); % y is correct
fz = @(z) S(z)*(y(z)/C(z))^(1.5)+  A*sqrt(y(z))-sqrt(mu)*dt; % f is correct

z0 = JoshBisection(fz,[-10,100]);
z = fzero(fz,z0);

flag = 1-(y(z)/norm(r1));
glag = A*sqrt(y(z)/mu);
% fdot = (sqrt(mu_e)/(norm(r1)*norm(r2)))*sqrt(y(z)/C(z))*(z*S(z)-1);
gdotlag = (1-(y(z)/norm(r2)));


end

