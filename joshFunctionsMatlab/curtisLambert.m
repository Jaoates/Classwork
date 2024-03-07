function [V1,V2,F,Fp,Fpz0,z,y,r1,r2] = curtisLambert(R1, R2, dt, string)
mu = 398600;

[C,S] = joshStumpffCoeffs();
C = @(z) sum(C.*joshStumpffZ(z));
S = @(z) sum(S.*joshStumpffZ(z));

%...Magnitudes of R1 and R2:
r1 = norm(R1);
r2 = norm(R2);
c12 = cross(R1, R2);
theta = acos(dot(R1,R2)/r1/r2);
%...Determine whether the orbit is prograde or retrograde:
if nargin < 4 || (~strcmp(string,'retro') & (~strcmp(string,'pro')))
    string = 'pro';
    fprintf('\n ** Prograde trajectory assumed.\n')
end
if strcmp(string,'pro')
    if c12(3) <= 0
        theta = 2*pi - theta;
    end
elseif strcmp(string,'retro')
    if c12(3) >= 0
        theta = 2*pi - theta;
    end
end
%...Equation 5.35:
A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));


disp("curt")
disp(theta)
disp(r1)
disp(r2)
disp(A)



%...Determine approximately where Finternal(z,t) changes sign, and
%...use that value of z as the starting value for Equation 5.45:
z = -100;
while Finternal(z,dt) < 0
    z = z + 0.1;
end
%...Set an error tolerance and a limit on the number of iterations:
tol = 1.e-8;
nmax = 5000;
%...Iterate on Equation 5.45 until z is determined to within the
%...error tolerance:
ratio = 1;
n = 0;
while (abs(ratio) > tol) & (n <= nmax)
    n = n + 1;
    ratio = Finternal(z,dt)/dFdz(z);
    z = z - ratio;
end
%...Report if the maximum number of iterations is exceeded:
if n >= nmax
    fprintf('\n\n **Number of iterations exceeds %g \n\n ',nmax)
end
%...Equation 5.46a:
f = 1 - yinternal(z)/r1;
%...Equation 5.46b:
g = A*sqrt(yinternal(z)/mu);
%...Equation 5.46d:
gdot = 1 - yinternal(z)/r2;
%...Equation 5.28:
V1 = 1/g*(R2 - f*R1);
%...Equation 5.29:
V2 = 1/g*(gdot*R2 - R1);


F = @(z)(yinternal(z)/C(z))^1.5*S(z) + A*sqrt(yinternal(z)) - sqrt(mu)*dt;
Fp = @(z)(yinternal(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...
                + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(yinternal(z)) ...
                + A*sqrt(C(z)/yinternal(z)));
Fpz0 = @(z)sqrt(2)/40*yinternal(0)^1.5 + A/8*(sqrt(yinternal(0)) + A*sqrt(1/2/yinternal(0)));
y = @(z) yinternal(z);




return
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Subfunctions used in the main body:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Equation 5.38:
    function dum = yinternal(z)
        dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
    end
%...Equation 5.40:
    function dum = Finternal(z,t)
        dum = (yinternal(z)/C(z))^1.5*S(z) + A*sqrt(yinternal(z)) - sqrt(mu)*t;
    end
%...Equation 5.43:
    function dum = dFdz(z)
        if z == 0
            dum = sqrt(2)/40*yinternal(0)^1.5 + A/8*(sqrt(yinternal(0)) + A*sqrt(1/2/yinternal(0)));
        else
            dum = (yinternal(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...
                + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(yinternal(z)) ...
                + A*sqrt(C(z)/yinternal(z)));
        end
    end
%...Stumpff functions:

end %lambert
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~