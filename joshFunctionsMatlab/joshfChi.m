function [fX,fpX] = joshfChi(r0,v0,mu,a,dt,C,S,pr)
% note fX = sqrt(mu)*dt - integral(rdX)
% fpX = -r

% it is highly recommended that C and S are passed in
arguments
    r0 double {mustBeReal}
    v0 double {mustBeReal}
    mu (1,1) double {mustBeReal}    
    a (1,1) double {mustBeReal}
    dt (1,1) double {mustBeReal}
    C (1,:) double {mustBeReal} = nan
    S (1,:) double {mustBeReal} = nan
    pr (1,1) double {mustBeReal} = nan
end

if isnan(C)|isnan(S)
    [C,S] = joshStumpffCoeffs();
end

[m1,n1] = size(r0);
[m2,n2] = size(v0);

if joshIsOnes([m1 n1 m2 n2])
    if isnan(pr)
        throw(MException("joshfChi:invalidInput","If r0 and v0 are given as vectors, pr should be the dot product of them"))        
    end
    
elseif(~joshIsOnes([m1 n1] == [m2 n2]))|~((n1==1&m1==3)|(n1==3&m1==1))
    throw(MException("joshfChi:invalidInput","r0 and v0 must be either 1x3 vectors or scalars"))
else
    pr = dot(r0,v0);
end


coefs = length(C);
if length(S)~=coefs
    throw(MException("joshfChi:invalidInput","S and C should be the same length"))
end

sm = sqrt(mu);
r0 = norm(r0);
% v0 = norm(v0);

% fX = @(X) (pr/sm)*X^2*C(z) + (1-a*r0)*X^3*S(z) + r0*X - sm*dt;
% fpX = @(X) (pr/sm)*X*(1-a*X^2*S(z)) + (1-a*r0)*X^2*C(z) + r0;
% (sum(C.*joshStumpffZ(X^2/a)))

% fX = @(X) (pr/sm)*X^2*(sum(C.*joshStumpffZ(X^2*a,coefs))) + (1-a*r0)*X^3*(sum(S.*joshStumpffZ(X^2*a,coefs))) + r0*X - sm*dt;
% fpX = @(X) (pr/sm)*X*(1-a*X^2*(sum(S.*joshStumpffZ(X^2*a,coefs)))) + (1-a*r0)*X^2*(sum(C.*joshStumpffZ(X^2*a,coefs))) + r0;

a = 1/a;
fX = @(X) (pr/sm)*X^2*sum(C.*joshStumpffZ(X^2*a,coefs)) + (1-a*r0)*X^3*sum(S.*joshStumpffZ(X^2*a,coefs)) + r0*X - sm*dt;
fpX = @(X) (pr/sm)*X*(1-a*X^2*sum(S.*joshStumpffZ(X^2*a,coefs)))+(1-a*r0)*X^2*sum(C.*joshStumpffZ(X^2*a,coefs))+r0;



% X^2*(sum(S.*joshStumpffZ(((X^2)/a)))) = X^2*C(z)
%
%
end