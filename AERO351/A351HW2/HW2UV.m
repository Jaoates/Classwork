clear all;
close all;

addpath 'C:\joshFunctionsMatlab'
% t = 0;
% 
% r_e = 6378; % km
% zp = 2000; % km
% rp=zp+r_e; % km
% V = [];
% R = [];
% Xv = [];
% equalOne = [];
% a = -a;
% [fX,fpX] = joshfChi(r0,v0,mu_e,a,dt,Cc,Sc);
% [X, count, xVector, errorVector, errorRatioVector] = joshNewtons(fX,fpX,X0,1e-14); 
% Tf = 2*60*60;

mu_e = 398600;

coefs = 15;
[Cc,Sc] = joshStumpffCoeffs(coefs);

dt = 60*60;
a0 = -19655;%km
a = 1/a0;

v0 = 10;
vr0 = 3.0752;
r0 = 10000;
pr=vr0*r0;

% z = @(X) (X^2/a0);
sm = sqrt(mu_e);

X0 = sm*dt*abs(a);   

% C = @(z) sum(Cc.*joshStumpffZ(z,coefs));
% S = @(z) sum(Sc.*joshStumpffZ(z,coefs));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fXwrk = @(X) (pr/sm)*X^2*C(X^2/a0) + (1-a*r0)*X^3*S(X^2/a0) + r0*X - sm*dt;
% fpXwrk = @(X) (pr/sm)*X*(1-a*X^2*S(X^2/a0))+(1-a*r0)*X^2*C(X^2/a0)+r0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fX = @(X) (pr/sm)*X^2*sum(Cc.*joshStumpffZ(X^2*a,coefs)) + (1-a*r0)*X^3*sum(Sc.*joshStumpffZ(X^2*a,coefs)) + r0*X - sm*dt;
% fpX = @(X) (pr/sm)*X*(1-a*X^2*sum(Sc.*joshStumpffZ(X^2*a,coefs)))+(1-a*r0)*X^2*sum(Cc.*joshStumpffZ(X^2*a,coefs))+r0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fX,fpX] = joshfChi(r0,v0,mu_e,a0,dt,Cc,Sc,pr);

fX(X0)
fpX(X0)



% [fX,fpX] = joshfChi(r0,v0,mu_e,a0,dt,Cc,Sc);

% [X] = joshNewtons(fX,fpX,X0,1e-14); 
[X] = joshNewtons(fX,fpX,X0,1e-14)

% Ratio(X) = fX(X)/fpX(X)

% f = 1-(X^2/norm(r0))*C(X^2/a);
% g = dt-((1/sqrt(mu_e))*X^3*S(X^2/a));
% r = f*r0+g*v0;
% f_dot = (sqrt(mu_e)/(norm(r)*norm(r0)))*X*((X^2/a)*S(X^2/a)-1);
% g_dot = 1-(X^2/norm(r))*C(X^2/a);
% v = f_dot*r0 + g_dot*v0;
% 
% f*g_dot-f_dot*g



% function s = stumpS(z)
% if z > 0
% s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
% elseif z < 0
% s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
% else
% s = 1/6;
% end
% end
% 
% function c = stumpC(z)
% if z > 0
% c = (1 - cos(sqrt(z)))/z;
% elseif z < 0
% c = (cosh(sqrt(-z)) - 1)/(-z);
% else
% c = 1/2;
% end
% end

