Zp = 800;
mu = 398600;
ecc = .2;
rE = 6378;


%Find:
%Vp
%h
%ra
%va
%a
%parameter
%Period
%at 8791:
%   theta
%   v
%   gamma

rp = rE+Zp;
ra = 1.5*rp;
a = (ra+rp)/2;
h = (a*mu*(1-ecc^2))^.5;

vaz =@ (theta)  (mu/h)*(1+ecc*cosd(theta));
v_p = vaz(0);
v_a = vaz(180);

P = ((2*pi)/(mu^.5))*(a^(3/2));

Parameter = (h^2)/mu;