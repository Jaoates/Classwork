function [Mw,Ohmdot] = wheelEOM(X,Mc,GeoProperties)
w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);
q = [eta,eps'];
r = X(11:13);
v = X(14:16);
v = v*1000; % m
Ohm = X(17:19);

Is = GeoProperties.Is;

Mw = Mc+cross(w,Is*Ohm);

Ohmdot = inv(Is)*Mc;
end