function [T] = FEEDBACK_ODE(X,params)

etac = params.etac;
kd=params.kd;
kp=params.kp;
epsc=params.epsc;


w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);

% qc = conj(qc)     conj qc
epsc = -epsc;

% qe = conj(qc) (*) q    calc qe
epse = etac*eps+eta*epsc + joshCross(epsc)*eps;
etae = (etac*eta) - epsc'*eps;


T = kp.*epse + kd.*w;

end

