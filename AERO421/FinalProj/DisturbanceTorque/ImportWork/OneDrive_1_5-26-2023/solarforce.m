function F_s = solarforce(s, A, n, p) % s is solar vec, A is 6 long vector of area

    F_s = zeros(3,1);

    for i = 1:length(n)
        Sw = wetarea(s, n(:,i), A(i));
        F_s = F_s + (-p*s*Sw);
    end
end
