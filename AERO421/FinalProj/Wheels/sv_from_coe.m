%% Curtis
function [r, v] = sv_from_coe(h,ecc,raan,inc,aop,theta,mu)
% h = coe(1);
% ecc = coe(2);
% raan = coe(3);
% inc = coe(4);
% aop = coe(5);
% theta = coe(6);

rp = (h^2/mu) * (1/(1 + ecc*cos(theta))) * (cos(theta)*[1;0;0] + sin(theta)*[0;1;0]);
vp = (mu/h) * (-sin(theta)*[1;0;0] + (ecc + cos(theta))*[0;1;0]);
R3_W = [cos(raan) sin(raan) 0; -sin(raan) cos(raan) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc)];
R3_w = [cos(aop) sin(aop) 0; -sin(aop) cos(aop) 0; 0 0 1];
Q_pX = (R3_w*R1_i*R3_W)';
r = Q_pX*rp;
v = Q_pX*vp;
r = r';
v = v';
end