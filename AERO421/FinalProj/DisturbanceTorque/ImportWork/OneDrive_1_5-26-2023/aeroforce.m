function F_a = aeroforce(v, A, n, rho_a) % s is solar vec, A is 6 long vector of area

    F_a = zeros(3,1);

    for i = 1:length(n)
        Sw = wetarea(v/norm(v), n(:,i), A(i));
        F_a = F_a + (-rho_a*norm(v)^2 * v/norm(v) *Sw);
    end
end
