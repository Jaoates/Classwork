function [V,D] = joshBasisFixOLD(V,D,tol)
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

warning('testing and argument validation is far from complete')


% normalizes colunm vectors
V(:,1) = V(:,1)/norm(V(:,1));
V(:,2) = V(:,2)/norm(V(:,2));
V(:,3) = V(:,3)/norm(V(:,3));


temp1 = cross(V(:,1),V(:,2));
temp2 = V(:,3) - temp1;
temp3 = -V(:,3) - temp1;
if ~joshIsOnes(~isAlways(temp2 > tol))&~joshIsOnes(~isAlways(temp3 > tol))
    throw(MException('joshBasisFix:invalivInput','input basis is not close enough to orthogonal'))
end
clear temp1 temp2 temp3


% we have found the right V,D pair if
% isAlways(lamtemp(1)>lamtemp(2))&isAlways(lamtemp(2)>lamtemp(3))
% and
% joshIsRHONB(V)

% direction combinations
dircombos = [
[1 1 1];
[1 1 -1];
[1 -1 1];
[1 -1 -1];
[-1 1 1];
[-1 1 -1];
[-1 -1 1];
[-1 -1 -1];   
];

% order combinations
ordcombos = [
[1 2 3]
[1 3 2]
[2 1 3]
[2 3 1]
[3 1 2]
[3 2 1]
];

% Vtemp = zeros(3);
lam = [D(1,1),D(2,2),D(3,3)];
% lamtemp = zeros(1,3);

for k = 1:length(ordcombos)
    % permutes the order of the basis
    Vtemp1(:,1) = V(:,ordcombos(k,1));
    Vtemp1(:,2) = V(:,ordcombos(k,2));
    Vtemp1(:,3) = V(:,ordcombos(k,3));

    lamtemp1(1) = lam(ordcombos(k,1));
    lamtemp1(2) = lam(ordcombos(k,2));
    lamtemp1(3) = lam(ordcombos(k,3));
    


    % check the different direction combos
    for i = 1:length(dircombos)
        % if dircombos is -1 we flip it
        for j = 1:3
            Vtemp2(:,j) = Vtemp1(:,j)*dircombos(i,j);
            lamtemp2(j) = lamtemp1(j)*dircombos(i,j);
        end

        temp1(k,i) = isAlways(lamtemp2(1)>=lamtemp2(2))&isAlways(lamtemp2(2)>=lamtemp2(3));
        temp2(k,i) = joshIsRHONB(Vtemp2,tol);
        if temp1(k,i) & temp2(k,i)
            V = Vtemp2;
            D = [[lamtemp2(1),0,0];[0,lamtemp2(2),0];[0,0,lamtemp2(3)]];
            return
        end
    end

end

throw(MException('joshBasisFix:unexpected','Unexpectedly, a basis meeting the criteria was not found'))

end

