
function someDerivatives = orbit(t,Y,mu_e) % y is current state, t is time (unused), mu is passed in
x(1) = Y(1);
x(2) = Y(2);
x(3) = Y(3);

dx(1) = Y(4);
dx(2) = Y(5);
dx(3) = Y(6);

r = norm(x);

ddx(1) = -mu_e*x(1)/r^3;
ddx(2) = -mu_e*x(2)/r^3;
ddx(3) = -mu_e*x(3)/r^3;

someDerivatives = [dx' ; ddx'];
end