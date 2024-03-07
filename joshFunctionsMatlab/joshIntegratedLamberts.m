function [v1,v2] = joshIntegratedLamberts(r1,r2,dt,grade,mu)
% this function will call joshfLamberts but requires less thinking to use.
% the downside is that for situations where r and v vectors aren't known,
% This function cannot be used and because this function will generate
% stumpf functions for you (which is convient and i shouldve done a long time ago)
% it is computationally less effcient.

arguments
    r1 (3,1) double {mustBeReal, mustBeNonNan}
    r2 (3,1) double {mustBeReal, mustBeNonNan}
%     magOrVec {mustBeMember(magOrVec,{'magnitude','vector'})} = 'magnitude'
    dt (1,1) double {mustBePositive}
    grade string {mustBeMember(grade,{'pro','retro'})} = 'pro'
    mu (1,1) double {mustBePositive} = 398600
end



% ready to run lamberts from obj2 to obj3
theta = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
if strcmp(grade,'retro')
    theta = 2*pi - theta;
end
% set up stumpffys
% coefs = 15;
% [Cc,Sc]=joshStumpffCoeffs(coefs);
% C = @(z) sum(Cc.*joshStumpffZ(z,coefs));
% S = @(z) sum(Sc.*joshStumpffZ(z,coefs));


[fz,y,A,z,flag,glag,gdotlag] = joshfLambert(norm(r1),norm(r2),dt,theta,mu);%,Cc,Sc);
% [v1rb, v2rb] =lambertsRB(r1,r2,dt,1);

v1 = (1/glag)*(r2-flag*r1);
v2 = (1/glag)*(gdotlag*r2-r1);


end

