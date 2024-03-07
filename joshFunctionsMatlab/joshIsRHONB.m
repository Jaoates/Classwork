function isRHONB = joshIsRHONB(V,tol)
% checks if the V matrices provided are RH-ONB's within tol


arguments
    V
    tol = 1e-5
end

% warning('testing and argument validation is far from complete')

% checks if the vectors are unit vectors
temp = isAlways(1-norm(V(:,1))<tol) | isAlways(1-norm(V(:,2))<tol) | isAlways(1-norm(V(:,2))<tol);
if ~temp
    isRHONB = false;
    warning('Not all vectors are unit vectors')
    return
end

% checks if the basis is right handed
temp = cross(V(:,1),V(:,2));
temp = V(:,3) - temp;

if ~joshIsOnes(~isAlways(temp > tol))
    isRHONB = false;
    warning('v1 x v2 ~= v3')
    return
end


isRHONB = true;
return
end
