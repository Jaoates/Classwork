function [r,v,t,X] = joshOrbitCoastOde(r0,v0,dt,mu)
% a simple ode propegator for two body orbits
arguments
    r0 (3,1) double {mustBeReal,mustBeNonNan}
    v0 (3,1) double {mustBeReal,mustBeNonNan}
    dt (1,1) double {mustBeReal,mustBeNonNan}
    mu (1,1) double {mustBeReal,mustBeNonNan}= 398600;
end
X0 = [r0;v0];
options = odeset('RelTol', 1e-8,'AbsTol',1e-13);
[t,X] = ode45(@orbitODEFun,[0,dt],X0,options,mu);
Xe = X(end,:);
r = Xe(1:3)';
v = Xe(4:6)';
end

function Xdot = orbitODEFun(t,X,mu)
r = X(1:3);
v = X(4:6);
vdot = (-mu/norm(r)^3)*r;
rdot = v;
Xdot = [rdot;vdot];
end