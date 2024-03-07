%% Joshua Oates
% I took a uv sovler from https://www.mathworks.com/matlabcentral/fileexchange/44789-solve-lambert-s-problem-in-two-body-dynamics
% I also used the matlab central code at https://www.mathworks.com/matlabcentral/fileexchange/26348-robust-solver-for-lambert-s-orbital-boundary-value-problem
clear all
close all
clc
addpath("C:\joshFunctionsMatlab")
mu = 398600.4415

%% setup
r1 = [15945.34,0,0]
r2 = [12214.83899,10249.46731,0]
dt = 76*60

%% UV
disp("UV")

[v1,v2,iter] = lambert_universal(mu, 0, r1, dt, r2, 0);
[a,ecc,theta,inc,raan,aop,h,T,E]=joshCOE(r1,v1,mu,"vector");
ecc = norm(ecc);
disp("a: "+string(a)+" km")
disp("ecc: "+string(ecc))
disp("theta "+string(theta)+" rad")
disp("inc: "+string(inc)+" rad")
disp("raan: "+string(raan)+" rad")
disp("aop: "+string(aop)+" rad")

disp("Velocity vector 1 in km/s")
disp(v1)
%% Izzo Gooding
disp("Izzo Gooding")

[v1,v2,extremal_distances,exitflag] = lambert(r1,r2,dt/3600/24,0,mu); %#coder

[a,ecc,theta,inc,raan,aop,h,T,E]=joshCOE(r1,v1,mu,"vector");
ecc = norm(ecc);
disp("a: "+string(a)+" km")
disp("ecc: "+string(ecc))
disp("theta "+string(theta)+" rad")
disp("inc: "+string(inc)+" rad")
disp("raan: "+string(raan)+" rad")
disp("aop: "+string(aop)+" rad")

disp("Velocity vector 1 in km/s")
disp(v1)
%% min energy
disp("min energy")

[v1,tm1,tm2,tmabs] = lamberts_min(r1,r2,1,mu);

[a,ecc,theta,inc,raan,aop,h,T,E]=joshCOE(r1,v1,mu,"vector");
ecc = norm(ecc);
disp("a: "+string(a)+" km")
disp("ecc: "+string(ecc))
disp("theta "+string(theta)+" rad")
disp("inc: "+string(inc)+" rad")
disp("raan: "+string(raan)+" rad")
disp("aop: "+string(aop)+" rad")

disp("Velocity vector 1 in km/s")
disp(v1)
%% gauss

disp("gauss")

v1 = joshGauss(r1,r2,dt,1);

[a,ecc,theta,inc,raan,aop,h,T,E]=joshCOE(r1,v1,mu,"vector");
ecc = norm(ecc);
disp("a: "+string(a)+" km")
disp("ecc: "+string(ecc))
disp("theta "+string(theta)+" rad")
disp("inc: "+string(inc)+" rad")
disp("raan: "+string(raan)+" rad")
disp("aop: "+string(aop)+" rad")

disp("Velocity vector 1 in km/s")
disp(v1)


%% gauss fun
function v1 = joshGauss(r1,r2,dt,tm)
    mu = 398600.4415;
    rn1 = norm(r1);
    rn2 = norm(r2);
    cosdt = dot(r1, r2)/(rn1 * rn2);
    dth = asin(tm*sqrt(1 - cosdt^2));
    l = (rn1+rn2)/(4*sqrt(rn1*rn2)*(cos(dth/2)))-.5;
    m = (mu*dt^2)/(2*sqrt(rn1*rn2)*(cos(dth/2)))^3;
    y = [0,1];
    while(abs(y(end)-y(end-1)) > 1e-10)
        x = m/y(end)^2-l;
        X = (4/3)*(1+(6/5)*x+(6/5)*(8/7)*x^2+(6/5)*(8/7)*(10/9)*x^3);
        y = [y,1+X*(l+x)];

    end
    y= y(end);
    dE = 4*asin(sqrt(x));
    a = ((sqrt(mu)*dt)/(2*y*sqrt(rn1*rn2)*sin(dE/2)*cos(dth/2)))^2;
    f = 1-a/norm(r1)*(1-cos(dE));
    g = dt-a^(1.5)/sqrt(mu)*(dE-sin(dE));
    v1 = (r2-f*r1)/g;

end

%% functions

function [a,ecc,theta,inc,raan,aop,h,T,E] = joshCOE(R,V,u,magOrVec)

%%%%%%%%%%%%%%%%%%
% Revamped to do rads and fit new naming convention % 
%%%%%%%%%%%%%%%%%%


%COESOATESJOSHUA Takes postion and velocity vector and returns COES, all in
%ECI frame of reffrenece, km and seconds as units and degrees
% a = semi major axis
% ecc = eccentricity
% i = inclination
% raan = right accention acending node 
% aop = argument of periapsis
% theta = true anomaly

% will return T in s as a period
% and E (sometimes epsilon) in km^2/s^2 as specific mechanical energy 
% h is agular momentum 

% magOrVec is a parameter that can be set for vector inputs to return
% vector h and ecc

% scalar entry should only be used if the spacecraft is at apoapse or
% periapse

arguments
    R {mustBeNumeric, mustBeReal}
    V  {mustBeNumeric, mustBeReal}
    u (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 3.986004418 * (10^5) %km^3/s^2
    magOrVec {mustBeMember(magOrVec,{'magnitude','vector'})} = 'magnitude'
end

[m1,n1] = size(R);
[m2,n2] = size(V);

if joshIsOnes([m1 n1 m2 n2])
    magOrVec = 'magnitude';
    R = [0 0 R];
    V = [V 0 0];
    warning("joshCOE will assume that R and V are normal if the inputs are scalar ie: the craft is in a circular orbit or is at periapse or apoapse")
elseif (~joshIsOnes([m1 n1] == [m2 n2]))|~((n1==1&m1==3)|(n1==3&m1==1))
    throw(MException("COEsOatesJoshua:invalidInput","R and V must be either 1x3 vectors or scalars"))
end

% uearth = 3.986004418(8)	x 10^14 m^3/s^2
ihat=[1,0,0];
khat=[0,0,1];

Rm = norm(R);
Vm = norm(V);

%calculate orbital constants
h = cross(R,V); %angular momentum vector
hm = norm(h);
E = ( ( Vm^2 ) / 2) - ( u / Rm ) ; %specfic mechanical energy
%calculate COEs

a = -u / (2 * E ); %semi major axis in km
T = 2*pi*sqrt((a^3)/u); %period in s

e = (1/u) * (((Vm^2)-(u/Rm) ) * R - (dot(R,V) * V)); %eccentricity vector
em = norm (e); %magnitude of e

inc = acos((dot(khat,h))/hm); %inclination 
n = cross(khat,h); %node vector
nm = norm(n); %magnitude n

%raan
raan = acos(dot(ihat,n)/nm);
if n(2) < 0  %checks the vector relative to j to see if angle is positive or negative
    raan = 2*pi - raan;
end

%aop
aop = acos(dot(n,e)/(nm*em));
if e(3) < 0
    aop = 2*pi - aop;
end

%theta
theta = acos(dot(e,R)/(em*Rm));
if(dot(R,V) < 0) % cehck flight path angle to see if it is postive or negative
    theta = 2*pi -theta;
end

if strcmp(magOrVec, 'magnitude') %for magnitude mode, vectors will not be returned
    h = hm;
    e = em;
    if joshIsOnes([m1 n1 m2 n2]) % for scalar inputs it is not possible to calculate these values
        theta = NaN;
        inc = NaN;
        raan = NaN;
        aop = NaN;
    end
end
ecc = e;

end


% this function is adapted from the work of Reno Brown, Credit to him
function [v1,tm1,tm2,tmabs] = lamberts_min(r1, r2, tm, mu)
    cosdt = dot(r1, r2)/(norm(r1) * norm(r2));
    dt = asind(tm*sqrt(1 - cosdt^2));

    c = sqrt(norm(r1)^2 + norm(r2)^2 - 2*norm(r1)*norm(r2)*cosd(dt));
    s = (norm(r1) + norm(r2) + c)/2;
    a_min = s/2;
    p_min = (1 - cosd( dt))*norm(r1)*norm(r2)/c;
    e_min = sqrt(1 - 2*p_min/s);
    alphae = pi;
    beta_e = 2*asin(sqrt((s - c)/s));
    tm1 = sqrt(a_min^3/mu).*(alphae + (sin(beta_e)));
    tm2 = sqrt(a_min^3/mu).*(alphae -(sin(beta_e)));
    tmabs = (1/3)*sqrt(2/mu) * (s^(1.5) - (s - c)^(1.5));
    v1 = (r2 - (1 - norm(r2)/p_min*(1 - cosd( dt)))*r1 ) * sqrt(mu*p_min)/(norm(r1)*norm(r2)*sind( dt));
    % v1, t_min_e, t_min_abs
end


% Shiva Iyer (25 November, 2013)
% Calculate the Stumpff functions C(x) and S(x)
function [C,S] = stumpff(Z)
    C = []; S = [];
    for (z = Z)
        if (z > 0)
            sx = sqrt(z);
            C(end+1,1) = (1-cos(sx))/z;
            S(end+1,1) = (sx-sin(sx))/sx^3;
        elseif (z < 0)
            sx = sqrt(-z);
            C(end+1,1) = (1-cosh(sx))/z;
            S(end+1,1) = (sinh(sx)-sx)/sx^3;
        else
            C(end+1,1) = 1/2;
            S(end+1,1) = 1/6;
        end
    end
end

% Shiva Iyer (28 November, 2013)
% Solve Lambert's problem using universal variables. <r1>,<r2> are position
% vectors around the body <mu> at times <t1>,<t2>. If <path> == 1,solve for
% the long path (delta_nu >= pi), else for the short path (delta_nu < pi).
% Return initial and final velocity vectors <v1>,<v2> and iterations <iter>
function [v1,v2,iter] = lambert_universal(mu, t1, r1, t2, r2, path)
    nr1 = norm(r1);
    nr2 = norm(r2);
    mutt = sqrt(mu)*(t2-t1);
    dnu = acos(dot(r1, r2)/(nr1*nr2));
    A = sqrt(nr1*nr2*(1+cos(dnu)));
    if (path == 1)
        A = -A;
    end
    z = 0;
    v1 = zeros(3,1);
    v2 = zeros(3,1);
    for (iter = 1:50)
        [C,S] = stumpff(z);
        y = abs(nr1+nr2-A*(1-z*S)/sqrt(C));
        x = sqrt(y/C);
        t = x^3*S+A*sqrt(y);
        if (abs(t-mutt) < 1E-6)
            f = 1-y/nr1;
            g = A*sqrt(y/mu);
            gd = 1-y/nr2;
            v1 = (r2-f*r1)/g;
            v2 = (gd*r2-r1)/g;
            return;
        end
        if (abs(z) > 1E-6)
            Cp = (1-z*S-2*C)/(2*z);
            Sp = (C-3*S)/(2*z);
            tp = x^3*(Sp-1.5*S*Cp/C)+0.125*A*(3*S*sqrt(y)/C+A/x);
        else
            tp = (sqrt(2)/40)*y^1.5+0.125*A*(sqrt(y)+A*sqrt(0.5/y));
        end
        z = z-(t-mutt)/tp;
    end
end


% LAMBERT            Lambert-targeter for ballistic flights
%                    (Izzo, and Lancaster, Blanchard & Gooding)
%
% Usage:
%    [V1, V2, extremal_distances, exitflag] = lambert(r1, r2, tf, m, GM_central)
%
% Dimensions:
%             r1, r2 ->  [1x3]
%             V1, V2 ->  [1x3]
% extremal_distances ->  [1x2]
%              tf, m ->  [1x1]
%         GM_central ->  [1x1]
%
% This function solves any Lambert problem *robustly*. It uses two separate
% solvers; the first one tried is a new and unpublished algorithm developed
% by Dr. D. Izzo from the European Space Agency [1]. This version is extremely
% fast, but especially for larger [m] it still fails quite frequently. In such
% cases, a MUCH more robust algorithm is started (the one by Lancaster &
% Blancard [2], with modifcations, initial values and other improvements by
% R.Gooding [3]), which is a lot slower partly because of its robustness.
%
% INPUT ARGUMENTS:
% ======================================================================
%    name        units    description
% ======================================================================
%   r1, r1       [km]     position vectors of the two terminal points.
%     tf        [days]    time of flight to solve for
%      m          [-]     specifies the number of complete orbits to complete
%                         (should be an integer)
% GM_central   [km3/s2]   std. grav. parameter (G�M = mu) of the central body
%
% OUTPUT ARGUMENTS:
% ======================================================================
%   name             units   description
% ======================================================================
%  V1, V2             [km/s]  terminal velocities at the end-points
%  extremal_distances  [km]   minimum(1) and maximum(2) distance of the
%                             spacecraft to the central body.
%  exitflag             [-]   Integer containing information on why the
%                             routine terminated. A value of +1 indicates
%                             success; a normal exit. A value of -1
%                             indicates that the given problem has no
%                             solution and cannot be solved. A value of -2
%                             indicates that both algorithms failed to find
%                             a solution. This should never occur since
%                             these problems are well-defined, and at the
%                             very least it can be determined that the
%                             problem has no solution. Nevertheless, it
%                             still occurs sometimes for accidental
%                             erroneous input, so it provides a basic
%                             mechanism to check any application using this
%                             algorithm.
%
% This routine can be compiled to increase its speed by a factor of about
% 10-15, which is certainly advisable when the complete application requires
% a great number of Lambert problems to be solved. The entire routine is
% written in embedded MATLAB, so it can be compiled with the emlmex()
% function (older MATLAB) or codegen() function (MATLAB 2011a and later).
%
% To do this using emlmex(), make sure MATLAB's current directory is equal
% to where this file is located. Then, copy-paste and execute the following
% commands to the command window:
%
%    example_input = {...
%         [0.0, 0.0, 0.0], ...% r1vec
%         [0.0, 0.0, 0.0], ...% r2vec
%          0.0, ...           % tf
%          0.0, ...           % m
%          0.0};              % muC
%    emlmex -eg example_input lambert.m
%
% This is of course assuming your compiler is configured correctly. See the
% documentation of emlmex() on how to do that.
%
% Using codegen(), the syntax is as follows:
%
%    example_input = {...
%         [0.0, 0.0, 0.0], ...% r1vec
%         [0.0, 0.0, 0.0], ...% r2vec
%          0.0, ...           % tf
%          0.0, ...           % m
%          0.0};              % muC
%    codegen lambert.m -args example_input
%
% Note that in newer MATLAB versions, the code analyzer will complain about
% the pragma "%#eml" after the main function's name, and possibly, issue
% subsequent warnings related to this issue. To get rid of this problem, simply
% replace the "%#eml" directive with "%#codegen".
%
%
%
% References:
%
% [1] Izzo, D. ESA Advanced Concepts team. Code used available in MGA.M, on
%     http://www.esa.int/gsp/ACT/inf/op/globopt.htm. Last retreived Nov, 2009.
% [2] Lancaster, E.R. and Blanchard, R.C. "A unified form of Lambert's theorem."
%     NASA technical note TN D-5368,1969.
% [3] Gooding, R.H. "A procedure for the solution of Lambert's orbital boundary-value
%     problem. Celestial Mechanics and Dynamical Astronomy, 48:145�165,1990.
%
% See also lambert_low_ExpoSins.


% Please report bugs and inquiries to:
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Licence    : 2-clause BSD (see License.txt)


% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5

% If you want to cite this work in an academic paper, please use
% the following template:
%
% Rody Oldenhuis, orcid.org/0000-0002-3162-3660. "Lambert" <version>,
% <date you last used it>. MATLAB Robust solver for Lambert's
% orbital-boundary value problem.
% https://nl.mathworks.com/matlabcentral/fileexchange/26348



% -----------------------------------------------------------------
% Izzo's version:
% Very fast, but not very robust for more complicated cases
% -----------------------------------------------------------------
function [V1,...
          V2, ...
          extremal_distances,...
          exitflag] = lambert(r1vec,...
                              r2vec,...
                              tf,...
                              m,...
                              muC) %#coder

% original documentation:
%{
 This routine implements a new algorithm that solves Lambert's problem. The
 algorithm has two major characteristics that makes it favorable to other
 existing ones.

 1) It describes the generic orbit solution of the boundary condition
 problem through the variable X=log(1+cos(alpha/2)). By doing so the
 graph of the time of flight become defined in the entire real axis and
 resembles a straight line. Convergence is granted within few iterations
 for all the possible geometries (except, of course, when the transfer
 angle is zero). When multiple revolutions are considered the variable is
 X=tan(cos(alpha/2)*pi/2).

 2) Once the orbit has been determined in the plane, this routine
 evaluates the velocity vectors at the two points in a way that is not
 singular for the transfer angle approaching to pi (Lagrange coefficient
 based methods are numerically not well suited for this purpose).

 As a result Lambert's problem is solved (with multiple revolutions
 being accounted for) with the same computational effort for all
 possible geometries. The case of near 180 transfers is also solved
 efficiently.

  We note here that even when the transfer angle is exactly equal to pi
 the algorithm does solve the problem in the plane (it finds X), but it
 is not able to evaluate the plane in which the orbit lies. A solution
 to this would be to provide the direction of the plane containing the
 transfer orbit from outside. This has not been implemented in this
 routine since such a direction would depend on which application the
 transfer is going to be used in.

 please report bugs to dario.izzo@esa.int
%}

% adjusted documentation:
%{
 By default, the short-way solution is computed. The long way solution
 may be requested by giving a negative value to the corresponding
 time-of-flight [tf].

 For problems with |m| > 0, there are generally two solutions. By
 default, the right branch solution will be returned. The left branch
 may be requested by giving a negative value to the corresponding
 number of complete revolutions [m].
%}

% Authors
% .-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.
% Name       : Dr. Dario Izzo
% E-mail     : dario.izzo@esa.int
% Affiliation: ESA / Advanced Concepts Team (ACT)

% Made more readible and optimized for speed by Rody P.S. Oldenhuis
% Code available in MGA.M on   http://www.esa.int/gsp/ACT/inf/op/globopt.htm

% last edited 12/Dec/2009

% ADJUSTED FOR EML-COMPILATION 24/Dec/2009

    % initial values
    tol = 1e-14;    bad = false;     days = 86400;

    % work with non-dimensional units
    r1 = sqrt(r1vec*r1vec.');  r1vec = r1vec/r1;
    V = sqrt(muC/r1);          r2vec = r2vec/r1;
    T = r1/V;                  tf    = tf*days/T; % also transform to seconds

    % relevant geometry parameters (non dimensional)
    mr2vec = sqrt(r2vec*r2vec.');
    % make 100% sure it's in (-1 <= dth <= +1)
    dth = acos( max(-1, min(1, (r1vec*r2vec.')/mr2vec)) );

    % decide whether to use the left or right branch (for multi-revolution
    % problems), and the long- or short way
    leftbranch = sign(m);   longway = sign(tf);
    m = abs(m);             tf = abs(tf);
    if (longway < 0), dth = 2*pi - dth; end

    % derived quantities
    c      = sqrt(1 + mr2vec^2 - 2*mr2vec*cos(dth)); % non-dimensional chord
    s      = (1 + mr2vec + c)/2;                     % non-dimensional semi-perimeter
    a_min  = s/2;                                    % minimum energy ellipse semi major axis
    Lambda = sqrt(mr2vec)*cos(dth/2)/s;              % lambda parameter (from BATTIN's book)
    crossprd = [r1vec(2)*r2vec(3) - r1vec(3)*r2vec(2),...
                r1vec(3)*r2vec(1) - r1vec(1)*r2vec(3),...% non-dimensional normal vectors
                r1vec(1)*r2vec(2) - r1vec(2)*r2vec(1)];
    mcr       = sqrt(crossprd*crossprd.');           % magnitues thereof
    nrmunit   = crossprd/mcr;                        % unit vector thereof

    % Initial values
    % ---------------------------------------------------------

    % ELMEX requires this variable to be declared OUTSIDE the IF-statement
    logt = log(tf); % avoid re-computing the same value

    % single revolution (1 solution)
    if (m == 0)

        % initial values
        inn1 = -0.5233;      % first initial guess
        inn2 = +0.5233;      % second initial guess
        x1   = log(1 + inn1);% transformed first initial guess
        x2   = log(1 + inn2);% transformed first second guess

        % multiple revolutions (0, 1 or 2 solutions)
        % the returned soltuion depends on the sign of [m]
    else
        % select initial values
        if (leftbranch < 0)
            inn1 = -0.5234; % first initial guess, left branch
            inn2 = -0.2234; % second initial guess, left branch
        else
            inn1 = +0.7234; % first initial guess, right branch
            inn2 = +0.5234; % second initial guess, right branch
        end
        x1 = tan(inn1*pi/2);% transformed first initial guess
        x2 = tan(inn2*pi/2);% transformed first second guess
    end

    % since (inn1, inn2) < 0, initial estimate is always ellipse
    xx   = [inn1, inn2];  aa = a_min./(1 - xx.^2);
    bbeta = longway * 2*asin(sqrt((s-c)/2./aa));
    % make 100.4% sure it's in (-1 <= xx <= +1)
    aalfa = 2*acos(  max(-1, min(1, xx)) );

    % evaluate the time of flight via Lagrange expression
    y12  = aa.*sqrt(aa).*((aalfa - sin(aalfa)) - (bbeta-sin(bbeta)) + 2*pi*m);

    % initial estimates for y
    if m == 0
        y1 = log(y12(1)) - logt;
        y2 = log(y12(2)) - logt;
    else
        y1 = y12(1) - tf;
        y2 = y12(2) - tf;
    end

    % Solve for x
    % ---------------------------------------------------------

    % Newton-Raphson iterations
    % NOTE - the number of iterations will go to infinity in case
    % m > 0  and there is no solution. Start the other routine in
    % that case
    err = inf;  iterations = 0; xnew = 0;
    while (err > tol)
        % increment number of iterations
        iterations = iterations + 1;
        % new x
        xnew = (x1*y2 - y1*x2) / (y2-y1);
        % copy-pasted code (for performance)
        if m == 0, x = exp(xnew) - 1; else x = atan(xnew)*2/pi; end
        a = a_min/(1 - x^2);
        if (x < 1) % ellipse
            beta = longway * 2*asin(sqrt((s-c)/2/a));
            % make 100.4% sure it's in (-1 <= xx <= +1)
            alfa = 2*acos( max(-1, min(1, x)) );
        else % hyperbola
            alfa = 2*acosh(x);
            beta = longway * 2*asinh(sqrt((s-c)/(-2*a)));
        end
        % evaluate the time of flight via Lagrange expression
        if (a > 0)
            tof = a*sqrt(a)*((alfa - sin(alfa)) - (beta-sin(beta)) + 2*pi*m);
        else
            tof = -a*sqrt(-a)*((sinh(alfa) - alfa) - (sinh(beta) - beta));
        end
        % new value of y
        if m ==0, ynew = log(tof) - logt; else ynew = tof - tf; end
        % save previous and current values for the next iterarion
        % (prevents getting stuck between two values)
        x1 = x2;  x2 = xnew;
        y1 = y2;  y2 = ynew;
        % update error
        err = abs(x1 - xnew);
        % escape clause
        if (iterations > 15), bad = true; break; end
    end

    % If the Newton-Raphson scheme failed, try to solve the problem
    % with the other Lambert targeter.
    if bad
        % NOTE: use the original, UN-normalized quantities
        [V1, V2, extremal_distances, exitflag] = ...
            lambert_LancasterBlanchard(r1vec*r1, r2vec*r1, longway*tf*T, leftbranch*m, muC);
        return
    end

    % convert converged value of x
    if m==0, x = exp(xnew) - 1; else x = atan(xnew)*2/pi; end

    %{
      The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
      now need the conic. As for transfer angles near to pi the Lagrange-
      coefficients technique goes singular (dg approaches a zero/zero that is
      numerically bad) we here use a different technique for those cases. When
      the transfer angle is exactly equal to pi, then the ih unit vector is not
      determined. The remaining equations, though, are still valid.
    %}

    % Solution for the semi-major axis
    a = a_min/(1-x^2);

    % Calculate psi
    if (x < 1) % ellipse
        beta = longway * 2*asin(sqrt((s-c)/2/a));
        % make 100.4% sure it's in (-1 <= xx <= +1)
        alfa = 2*acos( max(-1, min(1, x)) );
        psi  = (alfa-beta)/2;
        eta2 = 2*a*sin(psi)^2/s;
        eta  = sqrt(eta2);
    else       % hyperbola
        beta = longway * 2*asinh(sqrt((c-s)/2/a));
        alfa = 2*acosh(x);
        psi  = (alfa-beta)/2;
        eta2 = -2*a*sinh(psi)^2/s;
        eta  = sqrt(eta2);
    end

    % unit of the normalized normal vector
    ih = longway * nrmunit;

    % unit vector for normalized [r2vec]
    r2n = r2vec/mr2vec;

    % cross-products
    % don't use cross() (emlmex() would try to compile it, and this way it
    % also does not create any additional overhead)
    crsprd1 = [ih(2)*r1vec(3)-ih(3)*r1vec(2),...
               ih(3)*r1vec(1)-ih(1)*r1vec(3),...
               ih(1)*r1vec(2)-ih(2)*r1vec(1)];
    crsprd2 = [ih(2)*r2n(3)-ih(3)*r2n(2),...
               ih(3)*r2n(1)-ih(1)*r2n(3),...
               ih(1)*r2n(2)-ih(2)*r2n(1)];

    % radial and tangential directions for departure velocity
    Vr1 = 1/eta/sqrt(a_min) * (2*Lambda*a_min - Lambda - x*eta);
    Vt1 = sqrt(mr2vec/a_min/eta2 * sin(dth/2)^2);

    % radial and tangential directions for arrival velocity
    Vt2 = Vt1/mr2vec;
    Vr2 = (Vt1 - Vt2)/tan(dth/2) - Vr1;

    % terminal velocities
    V1 = (Vr1*r1vec + Vt1*crsprd1)*V;
    V2 = (Vr2*r2n + Vt2*crsprd2)*V;

    % exitflag
    exitflag = 1; % (success)

    % also compute minimum distance to central body
    % NOTE: use un-transformed vectors again!
    extremal_distances = ...
        minmax_distances(r1vec*r1, r1, r2vec*r1, mr2vec*r1, dth, a*r1, V1, V2, m, muC);

end

% -----------------------------------------------------------------
% Lancaster & Blanchard version, with improvements by Gooding
% Very reliable, moderately fast for both simple and complicated cases
% -----------------------------------------------------------------
function [V1,...
          V2,...
          extremal_distances,...
          exitflag] = lambert_LancasterBlanchard(r1vec,...
                                                 r2vec,...
                                                 tf,...
                                                 m,...
                                                 muC) %#coder
%{
LAMBERT_LANCASTERBLANCHARD       High-Thrust Lambert-targeter

lambert_LancasterBlanchard() uses the method developed by
Lancaster & Blancard, as described in their 1969 paper. Initial values,
and several details of the procedure, are provided by R.H. Gooding,
as described in his 1990 paper.
%}

% Please report bugs and inquiries to:
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Licence    : 2-clause BSD (see License.txt)


% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5

    % ADJUSTED FOR EML-COMPILATION 29/Sep/2009

    % manipulate input
    tol     = 1e-12;                            % optimum for numerical noise v.s. actual precision
    r1      = sqrt(r1vec*r1vec.');              % magnitude of r1vec
    r2      = sqrt(r2vec*r2vec.');              % magnitude of r2vec
    r1unit  = r1vec/r1;                         % unit vector of r1vec
    r2unit  = r2vec/r2;                         % unit vector of r2vec
    crsprod = cross(r1vec, r2vec, 2);           % cross product of r1vec and r2vec
    mcrsprd = sqrt(crsprod*crsprod.');          % magnitude of that cross product
    th1unit = cross(crsprod/mcrsprd, r1unit);   % unit vectors in the tangential-directions
    th2unit = cross(crsprod/mcrsprd, r2unit);
    % make 100.4% sure it's in (-1 <= x <= +1)
    dth = acos( max(-1, min(1, (r1vec*r2vec.')/r1/r2)) ); % turn angle

    % if the long way was selected, the turn-angle must be negative
    % to take care of the direction of final velocity
    longway = sign(tf); tf = abs(tf);
    if (longway < 0), dth = dth-2*pi; end

    % left-branch
    leftbranch = sign(m); m = abs(m);

    % define constants
    c  = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(dth));
    s  = (r1 + r2 + c) / 2;
    T  = sqrt(8*muC/s^3) * tf;
    q  = sqrt(r1*r2)/s * cos(dth/2);

    % general formulae for the initial values (Gooding)
    % -------------------------------------------------

    % some initial values
    T0  = LancasterBlanchard(0, q, m);
    Td  = T0 - T;
    phr = mod(2*atan2(1 - q^2, 2*q), 2*pi);

    % initial output is pessimistic
    V1 = NaN(1,3);    V2 = V1;    extremal_distances = [NaN, NaN];

    % single-revolution case
    if (m == 0)
        x01 = T0*Td/4/T;
        if (Td > 0)
            x0 = x01;
        else
            x01 = Td/(4 - Td);
            x02 = -sqrt( -Td/(T+T0/2) );
            W   = x01 + 1.7*sqrt(2 - phr/pi);
            if (W >= 0)
                x03 = x01;
            else
                x03 = x01 + (-W).^(1/16).*(x02 - x01);
            end
            lambda = 1 + x03*(1 + x01)/2 - 0.03*x03^2*sqrt(1 + x01);
            x0 = lambda*x03;
        end

        % this estimate might not give a solution
        if (x0 < -1), exitflag = -1; return; end

    % multi-revolution case
    else

        % determine minimum Tp(x)
        xMpi = 4/(3*pi*(2*m + 1));
        if (phr < pi)
            xM0 = xMpi*(phr/pi)^(1/8);
        elseif (phr > pi)
            xM0 = xMpi*(2 - (2 - phr/pi)^(1/8));
        % EMLMEX requires this one
        else
            xM0 = 0;
        end

        % use Halley's method
        xM = xM0;  Tp = inf;  iterations = 0;
        while abs(Tp) > tol
            % iterations
            iterations = iterations + 1;
            % compute first three derivatives
            [dummy, Tp, Tpp, Tppp] = LancasterBlanchard(xM, q, m);%#ok
            % new value of xM
            xMp = xM;
            xM  = xM - 2*Tp.*Tpp ./ (2*Tpp.^2 - Tp.*Tppp);
            % escape clause
            if mod(iterations, 7), xM = (xMp+xM)/2; end
            % the method might fail. Exit in that case
            if (iterations > 25), exitflag = -2; return; end
        end

        % xM should be elliptic (-1 < x < 1)
        % (this should be impossible to go wrong)
        if (xM < -1) || (xM > 1), exitflag = -1; return; end

        % corresponding time
        TM = LancasterBlanchard(xM, q, m);

        % T should lie above the minimum T
        if (TM > T), exitflag = -1; return; end

        % find two initial values for second solution (again with lambda-type patch)
        % --------------------------------------------------------------------------

        % some initial values
        TmTM  = T - TM;   T0mTM = T0 - TM;
        [dummy, Tp, Tpp] = LancasterBlanchard(xM, q, m);%#ok

        % first estimate (only if m > 0)
        if leftbranch > 0
            x   = sqrt( TmTM / (Tpp/2 + TmTM/(1-xM)^2) );
            W   = xM + x;
            W   = 4*W/(4 + TmTM) + (1 - W)^2;
            x0  = x*(1 - (1 + m + (dth - 1/2)) / ...
                (1 + 0.15*m)*x*(W/2 + 0.03*x*sqrt(W))) + xM;

            % first estimate might not be able to yield possible solution
            if (x0 > 1), exitflag = -1; return; end

        % second estimate (only if m > 0)
        else
            if (Td > 0)
                x0 = xM - sqrt(TM/(Tpp/2 - TmTM*(Tpp/2/T0mTM - 1/xM^2)));
            else
                x00 = Td / (4 - Td);
                W = x00 + 1.7*sqrt(2*(1 - phr));
                if (W >= 0)
                    x03 = x00;
                else
                    x03 = x00 - sqrt((-W)^(1/8))*(x00 + sqrt(-Td/(1.5*T0 - Td)));
                end
                W      = 4/(4 - Td);
                lambda = (1 + (1 + m + 0.24*(dth - 1/2)) / ...
                    (1 + 0.15*m)*x03*(W/2 - 0.03*x03*sqrt(W)));
                x0     = x03*lambda;
            end

            % estimate might not give solutions
            if (x0 < -1), exitflag = -1; return; end

        end
    end

    % find root of Lancaster & Blancard's function
    % --------------------------------------------

    % (Halley's method)
    x = x0; Tx = inf; iterations = 0;
    while abs(Tx) > tol
        % iterations
        iterations = iterations + 1;
        % compute function value, and first two derivatives
        [Tx, Tp, Tpp] = LancasterBlanchard(x, q, m);
        % find the root of the *difference* between the
        % function value [T_x] and the required time [T]
        Tx = Tx - T;
        % new value of x
        xp = x;
        x  = x - 2*Tx*Tp ./ (2*Tp^2 - Tx*Tpp);
        % escape clause
        if mod(iterations, 7), x = (xp+x)/2; end
        % Halley's method might fail
        if iterations > 25, exitflag = -2; return; end
    end

    % calculate terminal velocities
    % -----------------------------

    % constants required for this calculation
    gamma = sqrt(muC*s/2);
    if (c == 0)
        sigma = 1;
        rho   = 0;
        z     = abs(x);
    else
        sigma = 2*sqrt(r1*r2/(c^2)) * sin(dth/2);
        rho   = (r1 - r2)/c;
        z     = sqrt(1 + q^2*(x^2 - 1));
    end

    % radial component
    Vr1    = +gamma*((q*z - x) - rho*(q*z + x)) / r1;
    Vr1vec = Vr1*r1unit;
    Vr2    = -gamma*((q*z - x) + rho*(q*z + x)) / r2;
    Vr2vec = Vr2*r2unit;

    % tangential component
    Vtan1    = sigma * gamma * (z + q*x) / r1;
    Vtan1vec = Vtan1 * th1unit;
    Vtan2    = sigma * gamma * (z + q*x) / r2;
    Vtan2vec = Vtan2 * th2unit;

    % Cartesian velocity
    V1 = Vtan1vec + Vr1vec;
    V2 = Vtan2vec + Vr2vec;

    % exitflag
    exitflag = 1; % (success)

    % also determine minimum/maximum distance
    a = s/2/(1 - x^2); % semi-major axis
    extremal_distances = minmax_distances(r1vec, r1, r1vec, r2, dth, a, V1, V2, m, muC);

end

% Lancaster & Blanchard's function, and three derivatives thereof
function [T, Tp, Tpp, Tppp] = LancasterBlanchard(x, q, m)

    % protection against idiotic input
    if (x < -1) % impossible; negative eccentricity
        x = abs(x) - 2;
    elseif (x == -1) % impossible; offset x slightly
        x = x + eps;
    end

    % compute parameter E
    E  = x*x - 1;

    % T(x), T'(x), T''(x)
    if x == 1 % exactly parabolic; solutions known exactly
        % T(x)
        T = 4/3*(1-q^3);
        % T'(x)
        Tp = 4/5*(q^5 - 1);
        % T''(x)
        Tpp = Tp + 120/70*(1 - q^7);
        % T'''(x)
        Tppp = 3*(Tpp - Tp) + 2400/1080*(q^9 - 1);

    elseif abs(x-1) < 1e-2 % near-parabolic; compute with series
        % evaluate sigma
        [sig1, dsigdx1, d2sigdx21, d3sigdx31] = sigmax(-E);
        [sig2, dsigdx2, d2sigdx22, d3sigdx32] = sigmax(-E*q*q);
        % T(x)
        T = sig1 - q^3*sig2;
        % T'(x)
        Tp = 2*x*(q^5*dsigdx2 - dsigdx1);
        % T''(x)
        Tpp = Tp/x + 4*x^2*(d2sigdx21 - q^7*d2sigdx22);
        % T'''(x)
        Tppp = 3*(Tpp-Tp/x)/x + 8*x*x*(q^9*d3sigdx32 - d3sigdx31);

    else % all other cases
        % compute all substitution functions
        y  = sqrt(abs(E));
        z  = sqrt(1 + q^2*E);
        f  = y*(z - q*x);
        g  = x*z - q*E;

        % BUGFIX: (Simon Tardivel) this line is incorrect for E==0 and f+g==0
        % d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) );
        % it should be written out like so:
        if (E<0)
            d = atan2(f, g) + pi*m;
        elseif (E==0)
            d = 0;
        else
            d = log(max(0, f+g));
        end

        % T(x)
        T = 2*(x - q*z - d/y)/E;
        %  T'(x)
        Tp = (4 - 4*q^3*x/z - 3*x*T)/E;
        % T''(x)
        Tpp = (-4*q^3/z * (1 - q^2*x^2/z^2) - 3*T - 3*x*Tp)/E;
        % T'''(x)
        Tppp = (4*q^3/z^2*((1 - q^2*x^2/z^2) + 2*q^2*x/z^2*(z - x)) - 8*Tp - 7*x*Tpp)/E;

    end
end

% series approximation to T(x) and its derivatives
% (used for near-parabolic cases)
function [sig, dsigdx, d2sigdx2, d3sigdx3] = sigmax(y)

    % preload the factors [an]
    % (25 factors is more than enough for 16-digit accuracy)
    persistent an;
    if isempty(an)
        an = [
            4.000000000000000e-001;     2.142857142857143e-001;     4.629629629629630e-002
            6.628787878787879e-003;     7.211538461538461e-004;     6.365740740740740e-005
            4.741479925303455e-006;     3.059406328320802e-007;     1.742836409255060e-008
            8.892477331109578e-010;     4.110111531986532e-011;     1.736709384841458e-012
            6.759767240041426e-014;     2.439123386614026e-015;     8.203411614538007e-017
            2.583771576869575e-018;     7.652331327976716e-020;     2.138860629743989e-021
            5.659959451165552e-023;     1.422104833817366e-024;     3.401398483272306e-026
            7.762544304774155e-028;     1.693916882090479e-029;     3.541295006766860e-031
            7.105336187804402e-033];
    end

    % powers of y
    powers = y.^(1:25);

    % sigma itself
    sig = 4/3 + powers*an;

    % dsigma / dx (derivative)
    dsigdx = ( (1:25).*[1, powers(1:24)] ) * an;

    % d2sigma / dx2 (second derivative)
    d2sigdx2 = ( (1:25).*(0:24).*[1/y, 1, powers(1:23)] ) * an;

    % d3sigma / dx3 (third derivative)
    d3sigdx3 = ( (1:25).*(0:24).*(-1:23).*[1/y/y, 1/y, 1, powers(1:22)] ) * an;

end


% -----------------------------------------------------------------
% Helper functions
% -----------------------------------------------------------------

% compute minimum and maximum distances to the central body
function extremal_distances = minmax_distances(r1vec, r1,...
                                               r2vec, r2,...
                                               dth,...
                                               a,...
                                               V1, V2,...
                                               m,...
                                               muC)

    % default - minimum/maximum of r1,r2
    minimum_distance = min(r1,r2);
    maximum_distance = max(r1,r2);

    % was the longway used or not?
    longway = abs(dth) > pi;

    % eccentricity vector (use triple product identity)
    evec = ((V1*V1.')*r1vec - (V1*r1vec.')*V1)/muC - r1vec/r1;

    % eccentricity
    e = sqrt(evec*evec.');
    % apses
    pericenter = a*(1-e);
    apocenter  = inf;                    % parabolic/hyperbolic case
    if (e < 1), apocenter = a*(1+e); end % elliptic case

    % since we have the eccentricity vector, we know exactly where the
    % pericenter lies. Use this fact, and the given value of [dth], to
    % cross-check if the trajectory goes past it
    if (m > 0) % obvious case (always elliptical and both apses are traversed)
        minimum_distance = pericenter;
        maximum_distance = apocenter;
    else % less obvious case
        % compute theta1&2 ( use (AxB)-(CxD) = (C·B)(D·A) - (C·A)(B·D) ))
        pm1 = sign( r1*r1*(evec*V1.') - (r1vec*evec.')*(r1vec*V1.') );
        pm2 = sign( r2*r2*(evec*V2.') - (r2vec*evec.')*(r2vec*V2.') );
        % make 100.4% sure it's in (-1 <= theta12 <= +1)
        theta1 = pm1*acos( max(-1, min(1, (r1vec/r1)*(evec/e).')) );
        theta2 = pm2*acos( max(-1, min(1, (r2vec/r2)*(evec/e).')) );
        % points 1&2 are on opposite sides of the symmetry axis -- minimum
        % and maximum distance depends both on the value of [dth], and both
        % [theta1] and [theta2]
        if (theta1*theta2 < 0)
            % if |th1| + |th2| = turnangle, we know that the pericenter was
            % passed
            if abs(abs(theta1) + abs(theta2) - dth) < 5*eps(dth)
                minimum_distance = pericenter;
            % this condition can only be false for elliptic cases, and
            % when it is indeed false, we know that the orbit passed
            % apocenter
            else
                maximum_distance = apocenter;
            end
        % points 1&2 are on the same side of the symmetry axis. Only if the
        % long-way was used are the min. and max. distances different from
        % the min. and max. values of the radii (namely, equal to the apses)
        elseif longway
            minimum_distance = pericenter;
            if (e < 1), maximum_distance = apocenter; end
        end
    end

    % output argument
    extremal_distances = [minimum_distance, maximum_distance];

end

%% Improvement to Gauss -- Vallado Method
% inputs
function [r1, r2, r3, v2] =AERO557GaussExtensionVallado(r2,v2,rho1,rho2,rho3,f1,g1,f3,g3,rhohat1,rhohat2,rhohat3,mu,M,tau1,tau3,Rsite1,Rsite2,Rsite3)
%...Initialize the iterative improvement loop and set error tolerance:
rho1_old = rho1; rho2_old = rho2; rho3_old = rho3;
diff1 = 1; diff2 = 1; diff3 = 1;
n = 0;
nmax = 1000;
tol = 1.e-8;
%...Iterative improvement loop:
while ((diff1 > tol) && (diff2 > tol) && (diff3 > tol)) && (n < nmax)
n = n+1;
%...Compute quantities required by universal kepler's equation:
ro = norm(r2);
vo = norm(v2);
%vro = dot(v2,r2)/ro;
a = 2/ro - vo^2/mu;
%...Solve universal Kepler's equation at times tau1 and tau3 for
% universal anomalies x1 and x3:
[~,~,ff1,gg1,~,~,~] = AERO557kepler( r2,v2, tau1);
[~,~,ff3,gg3,~,~,~] = AERO557kepler( r2,v2, tau3);
%...Update the f and g functions at times tau1 and tau3 by
% averaging old and new:
f1 = (f1 + ff1)/2;
f3 = (f3 + ff3)/2;
g1 = (g1 + gg1)/2;
g3 = (g3 + gg3)/2;
%...
c1 = g3/(f1*g3 - f3*g1);
c3 = -g1/(f1*g3 - f3*g1);
%...
cmat(1,1)= -c1;
c2 = -1;
cmat(2,1)= -c2;
cmat(3,1)= -c3;
qmat = (M*cmat); %lir above
%...
r1 = Rsite1 + qmat(1,1)*rhohat1/c1;
r2 = Rsite2 + qmat(2,1)*rhohat2/c2;
r3 = Rsite3 + qmat(3,1)*rhohat3/c3;
%...
v2 = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1);
%...Calculate differences upon which to base convergence:
diff1 = abs(qmat(1,1) - rho1_old);
diff2 = abs(qmat(2,1) - rho2_old);
diff3 = abs(qmat(3,1) - rho3_old);
%...Update the slant ranges:
rho1_old = rho1; rho2_old = rho2; rho3_old = rho3;
end
%...End iterative improvement loop
%% Kepler
% ------------------------------------------------------------------------------
%
% function kepler
%
% this function solves keplers problem for orbit determination and returns a
% future geocentric equatorial (ijk) position and velocity vector. the
% solution uses universal variables.
%
% author : david vallado 719-573-2600 22 jun 2002
%
% revisions
% vallado - fix some mistakes 13 apr 2004
%
% inputs description range / units
% ro - ijk position vector - initial km
% vo - ijk velocity vector - initial km / s
% dtsec - length of time to propagate s
%
% outputs :
% r - ijk position vector km
% v - ijk velocity vector km / s
% error - error flag 'ok', ...
%
% locals :
% f - f expression
% g - g expression
% fdot - f dot expression
% gdot - g dot expression
% xold - old universal variable x
% xoldsqrd - xold squared
% xnew - new universal variable x
% xnewsqrd - xnew squared
% znew - new value of z
% c2new - c2(psi) function
% c3new - c3(psi) function
% dtsec - change in time s
% timenew - new time s
% rdotv - result of ro dot vo
% a - semi or axis km
% alpha - reciprocol 1/a
% sme - specific mech energy km2 / s2
% period - time period for satellite s
% s - variable for parabolic case
% w - variable for parabolic case
% h - angular momentum vector
% temp - temporary real*8 value
% i - index
%
% coupling :
% mag - magnitude of a vector
% findc2c3 - find c2 and c3 functions
%
% references :
% vallado 2004, 95-103, alg 8, ex 2-4
%
% [r,v] = kepler ( ro,vo, dtsec );
% ------------------------------------------------------------------------------
function [r,v,f,g,fdot,gdot,xnew] = AERO557kepler( ro,vo, dtseco)
%function [r,v,errork] = kepler ( ro,vo, dtseco, fid );
% ------------------------- implementation -----------------
% set constants and intermediate printouts
small = 1.0e-10;
infinite = 999999.9;
% ------------------------- mathematical --------------------
twopi = 2.0 * pi;
halfpi = pi * 0.5;
% ------------------------- conversions ---------------------
ft2m = 0.3048;
nm2m = 1852;
% ----------------------- physical constants ----------------
% WGS-84/EGM-96 constants used here
re = 6378.137; % km
flat = 1.0/298.257223563;
omegaearth = 7.292115e-5; % rad/s
mue = 398600.4418; % km3/s2
mum = 3.986004418e14; % m3/s2
% derived constants from the base values
eccearth = sqrt(2.0*flat - flat^2);
eccearthsqrd = eccearth^2;
renm = re / nm2m;
reft = re * 1000.0 / ft2m;
tusec = sqrt(re^3/mue);
tumin = tusec / 60.0;
tuday = tusec / 86400.0;
omegaearthradptu = omegaearth * tusec;
omegaearthradpmin = omegaearth * 60.0;
velkmps = sqrt(mue / re);
velftps = velkmps * 1000.0/ft2m;
velradpmin = velkmps * 60.0/re;
%for afspc
%velkmps1 = velradpmin*6378.135/60.0 7.90537051051763
%mu1 = velkmps*velkmps*6378.135 3.986003602567418e+005
degpsec = (180.0 / pi) / tusec;
radpday = 2.0 * pi * 1.002737909350795;
speedoflight = 2.99792458e8; % m/s
au = 149597870.0; % km
earth2moon = 384400.0; % km
moonradius = 1738.0; % km
sunradius = 696000.0; % km
masssun = 1.9891e30;
massearth = 5.9742e24;
massmoon = 7.3483e22;
show = 'n';
numiter = 50;
if show =='y'
fprintf(1,' ro %16.8f %16.8f %16.8f ER \n',ro/re );
fprintf(1,' vo %16.8f %16.8f %16.8f ER/TU \n',vo/velkmps );
end
% -------------------- initialize values -------------------
ktr = 0;
xold = 0.0;
znew = 0.0;
errork = ' ok';
dtsec = dtseco;
mulrev = 0;
if ( abs( dtseco ) > small )
magro = mag( ro );
magvo = mag( vo );
rdotv= dot( ro,vo );
% ------------- find sme, alpha, and a ------------------
sme= ( (magvo^2)*0.5 ) - ( mue /magro );
alpha= -sme*2.0/mue;
if ( abs( sme ) > small )
a= -mue / ( 2.0 *sme );
else
a= infinite;
end
if ( abs( alpha ) < small ) % parabola
alpha= 0.0;
end
if show =='y'
fprintf(1,' sme %16.8f a %16.8f alp %16.8f ER \n',sme/(mue/re),a/re, alpha*re );
fprintf(1,' sme %16.8f a %16.8f alp %16.8f km \n',sme, a,alpha );
fprintf(1,' ktr xn psi r xn+1dtn \n' );
end
% ------------ setup initial guess for x ---------------
% ----------------- circle and ellipse -------------------
if ( alpha >= small )
period= twopi * sqrt( abs(a)^3.0/mue );
% ------- next if needed for 2body multi-rev ----------
if ( abs( dtseco ) > abs( period ) )
mulrev = floor(dtseco/period);
end;
if ( abs(alpha-1.0 ) > small )
xold = sqrt(mue)*dtsec * alpha;
else
% - first guess can't be too close. ie a circle, r=a
xold = sqrt(mue) * dtsec * alpha * 0.97;
end
else
% -------------------- parabola ---------------------
if ( abs( alpha ) < small )
h = cross( ro,vo );
magh = mag(h);
p= magh*magh/mue;
s= 0.5 * (halfpi - datan( 3.0 *sqrt( mue / (p*p*p) )*dtsec ) );
w= atan( tan( s )^(1.0 /3.0 ) );
xold = sqrt(p) * ( 2.0 *cot(2.0 *w) );
alpha= 0.0;
else
% ------------------ hyperbola ------------------
temp= -2.0 * mue * dtsec / ...
( a*( rdotv + sign(dtsec)*sqrt(-mue*a)* ...
(1.0 -magro*alpha) ) );
xold= sign(dtsec) * sqrt(-a) *log(temp);
end
end
ktr= 1;
dtnew = -10.0;
while ((abs(dtnew/sqrt(mue) - dtsec) >= small) && (ktr < numiter))
xoldsqrd = xold*xold;
znew = xoldsqrd * alpha;
% ------------- find c2 and c3 functions --------------
[c2new,c3new] = findc2c3( znew );
% ------- use a newton iteration for new values -------
rval = xoldsqrd*c2new + rdotv/sqrt(mue) *xold*(1.0 -znew*c3new)
+ ...
magro*( 1.0 - znew*c2new );
dtnew= xoldsqrd*xold*c3new + rdotv/sqrt(mue)*xoldsqrd*c2new + ...
magro*xold*( 1.0 - znew*c3new );
% ------------- calculate new value for x -------------
xnew = xold + ( dtsec*sqrt(mue) - dtnew ) / rval;
if show =='y'
fprintf(1,'%3i %11.7f %11.7f %11.7f %11.7f %11.7f \n', ...
ktr,xold,znew,rval,xnew,dtnew);
fprintf(1,'%3i %11.7f %11.7f %11.7f %11.7f %11.7f \n', ...
ktr,xold/sqrt(re),znew,rval/re,xnew/sqrt(re),dtnew/sqrt(mue));
end
ktr = ktr + 1;
xold = xnew;
end
if ( ktr >= numiter )
errork= 'knotconv';
fprintf(1,'not converged in %2i iterations \n',numiter );
%v = zeros(1:3);
%r = zeros(1:3);
for i= 1 : 3
v(i)= 0.0;
r(i)= v(i);
end
else
% --- find position and velocity vectors at new time --
xnewsqrd = xnew*xnew;
f = 1.0 - ( xnewsqrd*c2new / magro );
g = dtsec - xnewsqrd*xnew*c3new/sqrt(mue);
%v = zeros(1:3);
%r = zeros(1:3);
for i= 1 : 3
r(i)= f*ro(i) + g*vo(i);
end
magr = mag( r );
gdot = 1.0 - ( xnewsqrd*c2new / magr );
fdot = ( sqrt(mue)*xnew / ( magro*magr ) ) * ( znew*c3new-1.0 );
for i= 1 : 3
v(i)= fdot*ro(i) + gdot*vo(i);
end
mag( v );
temp= f*gdot - fdot*g;
if ( abs(temp-1.0 ) > 0.00001 )
errork= 'fandg';
end
if show =='y'
fprintf(1,'f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n',f, g,fdot, gdot );
fprintf(1,'f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n',f,g/tusec, fdot*tusec, gdot );
fprintf(1,'r1 %16.8f %16.8f %16.8f ER \n',r/re );
fprintf(1,'v1 %16.8f %16.8f %16.8f ER/TU \n',v/velkmps );
end
end % if
else
% ----------- set vectors to incoming since 0 time --------
%v = zeros(1:3);
%r = zeros(1:3);
for i=1:3
r(i)= ro(i);
v(i)= vo(i);
end
end
end %kepler
% ------------------------------------------------------------------------------
%
% function findc2c3
%
% this function calculates the c2 and c3 functions for use in the universal
% variable calculation of z.
%
% author : david vallado 719-573-2600 27 may 2002
%
% revisions
% -
%
% inputs description range / units
% znew - z variable rad2
%
% outputs :
% c2new - c2 function value
% c3new - c3 function value
%
% locals :
% sqrtz - square root of znew
%
% coupling :
% sinh - hyperbolic sine
% cosh - hyperbolic cosine
%
% references :
% vallado 2001, 70-71, alg 1
%
% [c2new,c3new] = findc2c3 ( znew );
% ------------------------------------------------------------------------------
function [c2new,c3new] = findc2c3 ( znew );
small = 0.000001;
% ------------------------- implementation -----------------
if ( znew > small )
sqrtz = sqrt( znew );
c2new = (1.0 -cos( sqrtz )) / znew;
c3new = (sqrtz-sin( sqrtz )) / ( sqrtz^3 );
else
if ( znew < -small )
sqrtz = sqrt( -znew );
c2new = (1.0 -cosh( sqrtz )) / znew;
c3new = (sinh( sqrtz ) - sqrtz) / ( sqrtz^3 );
else
c2new = 0.5;
c3new = 1.0 /6.0;
end
end
end
% ------------------------------------------------------------------------------
%
% function mag
%
% this function finds the magnitude of a vector. the tolerance is set to
% 0.000001, thus the 1.0e-12 for the squared test of underflows.
%
% author : david vallado 719-573-2600 30 may 2002
%
% revisions
% vallado - fix tolerance to match coe, eq, etc 3 sep 2002
%
% inputs description range / units
% vec - vector
%
% outputs :
% mag - magnitude
%
% locals :
% none.
%
% coupling :
% none.
%
% mag = ( vec );
% ----------------------------------------------------------------------------- }
function mag = mag ( vec )
temp= vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3);
if abs( temp ) >= 1.0e-16
mag= sqrt( temp );
else
mag= 0.0;
end
end
end