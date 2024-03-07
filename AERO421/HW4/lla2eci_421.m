function pos = lla2eci(lat, lon, alt, JD)
% lat and long in degrees


Re = 6378.137;
gst = ct2lst(lon, JD); %deg
theta = (gst) * pi/180.0; %rad

r = (Re + alt)*cos(lat*pi/180.0); % km

pos(1) = r*cos(theta);% km
pos(2) = r*sin(theta);% km
pos(3) = (alt+Re)*sin(lat*pi/180.0);% km