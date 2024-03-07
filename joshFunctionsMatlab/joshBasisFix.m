function [V,D] = joshBasisFix(V,D,tol)
% takes a symbolic or double V and D as given by eig()
% returns a pair V and D such that the column vectors in V are eiganvectors
% and are a right handed orthonormal basis and the corresponding
% eiganvalues in the diagonal of D are in decending order.
% tol may have to be set to 0 for correct symbolic comparisons.


arguments
    V
    D
    tol = 1e-10
end


% possible permutations of order of D and consequently V
ordcombo = [
[1 2 3]
[1 3 2]
[2 1 3]
[2 3 1]
[3 1 2]
[3 2 1]
];

% for internal use
lam = [D(1,1),D(2,2),D(3,3)];

% permute until in order
for i = 1:length(ordcombo)
    % reorder lam
    lamtemp = [lam(ordcombo(i,1)),lam(ordcombo(i,2)),lam(ordcombo(i,3))];
    if isAlways(lamtemp(1)>=lamtemp(2))&isAlways(lamtemp(2)>=lamtemp(3))
        lam = lamtemp;
        % reorder V
        V = [V(:,ordcombo(i,1)),V(:,ordcombo(i,2)),V(:,ordcombo(i,3))];
        % break if correct order is found
        break
    end
end

% if we fail to prove order, throw error
if ~(isAlways(lam(1)>=lam(2))&isAlways(lam(2)>=lam(3)))
    throw(MException('joshBasisFix:Unexpected','Unexpectedly, no order of D could be found where the first value is greater than the second and the second is greater than the third'))
end

% normalizes colunm vectors
V(:,1) = V(:,1)/norm(V(:,1));
V(:,2) = V(:,2)/norm(V(:,2));
V(:,3) = V(:,3)/norm(V(:,3));

if ~joshIsRHONB(V,tol)
    % repair by flipping one value
    V(:,3) = -V(:,3);
    % check repair worked
    if ~joshIsRHONB(V,tol)
        throw(MException('joshBasisFix:Unexpected','Unexpectedly, the basis V could not be made into a right handed ONB'))
    end
end


% put back in form of D
D = [[lam(1),0,0];[0,lam(2),0];[0,0,lam(3)]];

return