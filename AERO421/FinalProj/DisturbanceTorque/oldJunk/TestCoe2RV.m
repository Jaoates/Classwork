% test script for sv to coe

mu = 398600;
h=53335.2
ecc=0
raan = 0
inc = deg2rad(98.43)
aop = 0
theta = 0
a = (h^2/mu)/(1-ecc)


% [r,v]  = keplerian2ijk(a,ecc,inc,raan,aop,theta)
[r,v]  = joshCOE2rv(a,ecc,theta,inc,raan,aop,mu)
[r,v] = sv_from_coe(h,ecc,raan,inc,aop,theta,mu)




