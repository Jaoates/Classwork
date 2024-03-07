function [isRotM] = joshIsRotM(M,tol)
    arguments
    M
    tol = 1e-14
    end
    decimals = -log10(tol);
    isRotM = (round(M*M',decimals) == eye(3) & round(M'*M,decimals) == eye(3) & round(det(M),decimals) == 1);
end

