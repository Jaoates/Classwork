function Xdot = EOM_ODElinear(t,X,I,kp,kd,epsc,etac)

eps = X(1:3);
epsdot = X(4:6);
% eta = sqrt(1-eps'*eps);
eta = 1;

% qc = conj(qc)     conj qc
% epsc = -epsc;

% qe = conj(qc) (*) q    calc qe

% epse = etac*eps + eta*epsc + joshCross(epsc)*eps;
% etae = (etac*eta) - epsc'*eps;

qc = [etac,epsc'];
q = [eta,eps'];

qe = quatmultiply(quatconj(qc),q);

etae = qe(1);
epse = qe(2:4)';

w = 2*epsdot;
% epsdotdot = inv(I)*(kp.*epse + kd.*2.*eps);
epsdotdot = inv(I)*(kp.*epse + kd.*2.*epsdot);

Xdot = zeros(6,1);
Xdot = [epsdot;epsdotdot];
end