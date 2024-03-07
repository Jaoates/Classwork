function Ta = AeroTorque(X,GeoProperties)

rho_a = 1.0135e-14;

w = X(1:3);
E = X(4:6);
eps = X(7:9);
eta = X(10);
q = [eta,eps'];
r = X(11:13);
v = X(14:16);
v = v*1000; % m

A = GeoProperties.A;
rho = GeoProperties.rho;

v_b = quatrotate(quatconj(q),v')';

n = [[1; 0; 0 ], [0; 1; 0], [0; 0; 1], [-1; 0; 0 ], [0; -1; 0], [0; 0; -1]]; %normal vecs
Cp = wettedcenter(v_b, A, n, rho);
F_a = aeroforce(v_b, A, n, rho_a);
Ta = cross(Cp, F_a);


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

function F_a = aeroforce(v, A, n, rho_a) % s is solar vec, A is 6 long vector of area

    F_a = zeros(3,1);

    for i = 1:length(n)
        Sw = wetarea(v/norm(v), n(:,i), A(i));
        F_a = F_a + (-rho_a*norm(v)^2 * v/norm(v) *Sw);
    end
end