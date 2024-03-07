function Xdot = EOM_ODE(t,X,I,T)
w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);

r = X(11:13);
v = X(14:16);

mu = 398600.4418;

% this matrix is part of the EOM for a 321 sequence for Edot. It is not a
% rotation matrix
mat = [[cos(E(2)),      sin(E(2))*sin(E(1)),            sin(E(2))*cos(E(1))];...
    [0,                 cos(E(2))*cos(E(1)),           -cos(E(2))*sin(E(1))];...
    [0,                 sin(E(1)),                      cos(E(1))           ]]*(1/cos(E(2)));


Edot = mat*w;

epsx = joshCross(eps);
epsdot = .5*(eta*eye(3)+epsx)*w;
etadot = -.5*eps'*w;

wx = joshCross(w);
wdot = -inv(I)*(wx*I*w-T);

rdot = v;
vdot = -mu*r/norm(r)^3;

Xdot = [wdot;Edot;epsdot;etadot;rdot;vdot];
end