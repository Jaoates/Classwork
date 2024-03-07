function Sw = wetarea(s, n, A)

    ndot = dot(n, s); % check if 'wetted'
        if ndot < 0
            ndot = 0; 
        end
    Sw = ndot * A;
end