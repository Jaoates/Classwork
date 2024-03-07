% function [v1,v2,z,f,g,fdot,gdot]  = joshLamberts(r1,r2,dt,mu,grade)
% arguments
%     r1 double {mustBeReal}
%     r2 double {mustBeReal}
%     dt (1,1) double {mustBePositive}
%     mu (1,1) double {mustBeReal} = 398600
%     grade {mustBeMember(grade,{'pro','retro'})} = 'pro'
% 
% end

% warning("joshLambert: This function may be useful but it is not well tested and complete argument validation has not been implimented.")
clear
r1 = [5000 10000 2100]
r2= [-14600 2500 7000]
grade = 'retro'

% r1 = [5644,2830,4170]; % km
% r2 = [-2240,7320,4980]; % km
% dt = 20; % min
% dt = dt*60; % sec

coefs = 15;
[Cc,Sc]=joshStumpffCoeffs(coefs);
C = @(z) sum(Cc.*joshStumpffZ(z,coefs));
S = @(z) sum(Sc.*joshStumpffZ(z,coefs));

mu_e = 398600;
dt = 3600; % sec
zc = cross(r1,r2);
zc = zc(3);

theta1 = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
theta2 = 2*pi - acos(dot(r1,r2)/(norm(r1)*norm(r2)));

% if (zc<0 & strcmp(grade,"pro")) | (zc>=0 & strcmp(grade,"retro"))
%     theta = 2*pi-theta;
% end


y = @(z)(r1+r2+  A*((z*S(z)-1))/sqrt(C(z)));

% [f,z,y,A] = joshfLambert(r1,r2,dt,theta,mu,C,S)
[f,zjosh,y,A] = joshfLambert(norm(r1),norm(r2),dt,theta2,mu_e,Cc,Sc);
[V1,V2,F,Fp,Fpz0,zcurtis,Y,r1curt,r2curt] = curtisLambert(r1,r2,dt,grade);


% z0 = JoshBracket(F,0,.00001,1,[-100,100]);
% z0 = mean(z0);
% 
% fp_Newtons = @(z,fp,fpz0)fp_zcheck(z,fp,fpz0);
% [z] = joshNewtons(f,fp_Newtons(z,fp,fpz0),z0,1e-14,200000)

z = fzero(f,0)

ffun = 1-(y(z)/norm(r1));
gfun = A*sqrt(y(z)/mu_e);
% fdot = (sqrt(mu_e)/(norm(r1)*norm(r2)))*sqrt(y(z)/C(z))*(z*S(z)-1);
gdot = (1-(y(z)/norm(r2)));

r1 = [5000 10000 2100]
r2=[-14600 2500 7000]

v1 = (1/gfun)*(r2-ffun*r1)
v2 = (1/gfun)*(gdot*r2-r1)
V1
V2


% functions
% for lamberts
function out = fp_zcheck(z,fp,fpz0)
if z == 0
    out = fpz0;
else
    out = fp;
end
end

% end


