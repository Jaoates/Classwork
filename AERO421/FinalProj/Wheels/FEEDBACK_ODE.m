function [T] = FEEDBACK_ODE(X,params)

% etac = params.etac;
kd=params.kd;
kp=params.kp;
% epsc=params.epsc;

%{
qbE ~ CbE
qLE ~ CLE
rb = qbE * rE * conj(qbE)
qbL = f(qbE,qLE)
rb = qbL * rL * conj(qbL)
rL = qLE * rE * conj(qLE)
rE = conj(qbE) * rb * qbE
rb = qbL * ( qLE * (conj(qbE) * rb * qbE) * conj(qLE) ) * conj(qbL)
rL = qLE * (conj(qbE) * rb * qbE) * conj(qLE)
rL = qLb * rb *conj(qLb)
qLb = qLE * conj(qbE)
qbL = conj(qLb)
%}

w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);
qbE = [eta,eps'];

r = X(11:13);
v = X(14:16);

qLE = RV2LVLH(r,v);
% qbL = quatconj(quatmultiply(qLE,quatconj(qbE)))


qc = qLE;
q = qbE;
qe = quatmultiply(quatconj(qc), q);

epse = qe(2:4)'; % the distance from b to L is the error

wLE = cross(r,v)/norm(r)^2;
wbE = w;
wbL = wbE-quatrotate(qbE,wLE')';

T = kp.*epse + kd.*wbL;
% T = kp.*epse + kd.*w;

end

