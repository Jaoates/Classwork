function [T_s] = SolarPressureTorque(X,GeoProperties)

p = 4.5e-6; %N/m^2
s_eci = [1; 0; 0];

w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);
q = [eta,eps'];
r = X(11:13);
v = X(14:16);

A = GeoProperties.A;
rho = GeoProperties.rho;

s_b = quatrotate(quatconj(q),s_eci')';

n = [[1; 0; 0 ], [0; 1; 0], [0; 0; 1], [-1; 0; 0 ], [0; -1; 0], [0; 0; -1]]; %normal vecs
Cp = wettedcenter(s_b, A, n, rho);
Fs = solarforce(s_b, A, n, p);
T_s = cross(Cp, Fs);

end

function C_ps = wettedcenter(s, A, n, rho) % for cubeish spacecraft aligned with body frame 

    C_ps = zeros(3,1);
    assert(length(A) == length(n));
    assert(length(rho) == length(n));

    for i = 1:length(n)
        Sw = wetarea(s, n(:,i), A(i));
        if Sw ~=0
           C_ps = C_ps + (rho(:,i) * Sw / Sw);
        end
    end
end

function Sw = wetarea(s, n, A)

    ndot = dot(n, s); % check if 'wetted'
        if ndot < 0
            ndot = 0; 
        end
    Sw = ndot * A;
end

function F_s = solarforce(s, A, n, p) % s is solar vec, A is 6 long vector of area

    F_s = zeros(3,1);

    for i = 1:length(n)
        Sw = wetarea(s, n(:,i), A(i));
        F_s = F_s + (-p*s*Sw);
    end
end
