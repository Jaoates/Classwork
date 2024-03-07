function [v1,v2,z,f,g,fdot,gdot]  = joshLamberts(r1,r2,dt,mu,grade)
% this function is defunct
arguments
    r1 double {mustBeReal}
    r2 double {mustBeReal}
    dt (1,1) double {mustBePositive}
    mu (1,1) double {mustBeReal} = 398600
    grade {mustBeMember(grade,{'pro','retro'})} = 'pro'

end
% function for refrence only
throw(MException("joshLambert:notSupported","This function is not supported and shouldn't be used. Use as a refrence only."))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning("joshLambert: This function may be useful but it is not well tested and complete argument validation has not been implimented.")

coefs = 15;
[Cc,Sc]=joshStumpffCoeffs(coefs);
C = @(z) sum(Cc.*joshStumpffZ(z,coefs));
S = @(z) sum(Sc.*joshStumpffZ(z,coefs));

mu_e = 398600;
dt = 3600; % sec
zc = cross(r1,r2);
zc = zc(3);

theta = acos(dot(r1,r2)/(norm(r1)*norm(r2)));

if (zc<0 & strcmp(grade,"pro")) | (zc>=0 & strcmp(grade,"retro"))
    theta = 2*pi-theta;
end

clear theta2 theta1

y = @(z)(r1+r2+  A*((z*S(z)-1))/sqrt(C(z)));

[f,fp,fpz0] = joshfLambert(norm(r1),norm(r2),mu_e,dt,theta,Cc,Sc);

z0 = JoshBracket(f,0.1,.001,1,[-2,2]);
z0 = mean(z0);

fp_Newtons = @(z,fp,fpz0)fp_zcheck(z,fp,fpz0);
[z] = joshNewtons(f,fp_Newtons(z,fp,fpz0),z0,1e-14); 



% functions
% for lamberts
function out = fp_zcheck(z,fp,fpz0)
if z == 0
    out = fpz0;
else
    out = fp;
end
end

end


