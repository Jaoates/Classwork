function [v1, v2] = lambertsRB(r1, r2, dt, pro)

coefs = 15;
[Cc,Sc]=joshStumpffCoeffs(coefs);
C = @(z) sum(Cc.*joshStumpffZ(z,coefs));
S = @(z) sum(Sc.*joshStumpffZ(z,coefs));

    muearth = 398600;
    dtheta = acosd(dot(r1,r2) / (norm(r1) * norm(r2)));
    test = cross(r1, r2);
    if pro * test(3) < 0
        dtheta = 360 -dtheta;
    end

    A = sind(dtheta) * sqrt(norm(r1) * norm(r2) / (1 - cosd(dtheta)));

    y = @(z) norm(r1) + norm(r2) + A*(z*S(z) - 1) / sqrt(C(z));
    F = @(z) (y(z)/C(z))^1.5 * S(z) + A * sqrt(y(z)) - sqrt(muearth)*dt;

    z = fzero(F, 0);

    f = 1 - y(z)/norm(r1);
    g = A * sqrt(y(z)/muearth);
    gdot = 1 - y(z) / norm(r2);

    v1 = 1/g * (r2 - f*r1);
    v2 = 1/g * (gdot*r2 - r1);
end