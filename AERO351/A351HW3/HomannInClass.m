za1 = 300;
zp1 = 300;
mu_e = 398600;
r_e = 6378; % km
r1 = za1+r_e;
% r1 = za2+r_e;
r2 = 3000+r_e;

ecc = (r2-r1)/(r2+r1);
v1 = sqrt(mu_e/r1);
v2 = sqrt(mu_e/r2);

h1 = r1*v1;
h2 = r2*v2;

ht = sqrt((r1*(1+ecc))*mu_e);

vtp = ht/r1;
vta = ht/r2;


