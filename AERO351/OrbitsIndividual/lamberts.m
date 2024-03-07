%% 5.6
function [v1,v2,z,f,g,fdot,gdot]  = joshLamberts(r1,r2,dt,mu,)
coefs = 15;
[Cc,Sc]=joshStumpffCoeffs(coefs);
C = @(z) sum(Cc.*joshStumpffZ(z,coefs));
S = @(z) sum(Sc.*joshStumpffZ(z,coefs));

mu_e = 398600;
r_e = 6378; % km

r1 = [5000 10000 2100]; % km
r2 =  [-14600 2500 7000]; % km

dt = 3600; % sec

z = cross(r1,r2);
z = z(3);

theta1 = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
theta2 = 2*pi - acos(dot(r1,r2)/(norm(r1)*norm(r2)));

y = @(z)(r1+r2+  A*((z*S(z)-1))/sqrt(C(z)));

[f,fp,fpz0] = joshfLambert(norm(r1),norm(r2),mu_e,dt,theta1,Cc,Sc);



z0 = JoshBracket(f,0.1,.001,1,[-2,2]);
z0 = mean(z0);
[z] = joshNewtons(f,fp_zcheck(z,fp,fpz0),z0,1e-14); 





% functions
% for lamberts
function out = fp_zcheck(z,fp,fpz0)
if z == 0
    out = fpz0;
else
    out = fp;
end
end




