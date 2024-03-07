function Tg  = gravgrad(X,I)

% break out X
w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);

r = X(11:13);
v = X(14:16);

R = r*1000;%m

% work out Tg
quat = quaternion([eta eps']);
C21 = rotmat(quat,'point');
C12 = C21';

rb = C12*R;

mu = 398600;%km
mu = 3.986004418e14;%m
Tg = ((3*mu)/(norm(rb)^5))*joshCross(rb)*I*rb;
end