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

