function [jd,thetaTime ,j2000, j0, ut, thetaG, jcent2000] = joshJulian(t,thetaLongitude)
% takes t as a datetime object in UT and a longitude
% jd - juliandate day
% thetaTime - local sidereal time in degrees
% j2000 - julian date from 2000 (jd-j2000_0)
% j0 - julian days
% ut - UT in hours
% thetaG - Grennich sidereal time
% number of julian centuries from j2000
arguments
    t (1,1) datetime
    thetaLongitude (1,1) {mustBeReal} = 0
end


j2000_0 = 2451545;

[yr,mo,da] = ymd(t);
[hr,mn,sc] = hms(t);

j0 = 367*yr-floor((7*(yr+floor((mo+9)/12)))/4)+floor((275*mo)/9)+da+1721013.5;
ut = hr + mn/60 + sc/3600;
jd = j0 + (ut/24);
j2000 = jd - j2000_0;

t0 = (j0 - j2000_0)/36525;

thetaG0 = 100.4606184 + 36000.77004*t0 + 0.000387933*t0^2 - 2.583*(10e-8)*t0^3;
thetaG = thetaG0 + 360.98564724* (ut/24);
thetaTime = thetaG + thetaLongitude;

thetaG = mod(thetaG,360);
thetaTime = mod(thetaTime,360);

jcent2000 = j2000/36525;

end
