function [dstate] = BIGorbitEQ(time, state, muearth)
%The core orbital eq broken into components

x=state(1);
y=state(2);
z=state(3);

dx=state(4);
dy=state(5);
dz=state(6);

r=norm([x y z]);

ddx=-muearth*x/r^3;
ddy=-muearth*y/r^3;
ddz=-muearth*z/r^3;

dstate=[dx; dy; dz; ddx; ddy; ddz];

end