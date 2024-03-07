function Tm = MagTorque(t,X,GeoProperties)

w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);
q = [eta,eps'];
r = X(11:13);
v = X(14:16);
v = v*1000; % m

r = r*1000;
m = GeoProperties.magneticDipole;

jd_0 = 2459946; % equinox?
jd = jd_0+t/86400;
UTC = datevec(datetime(jd, "ConvertFrom","juliandate"));
lla = eci2lla(r',UTC);
% [XYZ, H, D, I, F] = wrldmagm_fixed(height, lat, lon, dyear, varargin);
B = wrldmagm(lla(3),lla(1),lla(2),decyear(datetime(jd, "ConvertFrom","juliandate")));
B = quatrotate(q,B')'; % ECI to Body

B = B/1e9;

Tm = cross(m,B);


end