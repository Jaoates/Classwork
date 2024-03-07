function [PI,a] = joshBuckPiTheory(M,qi)
%JOSHBUCKPITHEORY finds pis of dimentional matrix M

% to find physically meaningful function h in terms of non-dimentional
% params pi1,pi2,...,pip
% h(q1,q2,...qn)=0
% H(pi1,pi2,...pip)=0

arguments
    M (:,:) double
    qi (:,1) string = []
end


a = null(M,'rational'); % a's, should match length q
[k,n] = size(M);

p = n-k; % number of pi's

if p<1
    throw(MException("joshBuckPiTheory:invalidInput","M is wrong dimensions."))  
end

for i = 1:p
    acol = a(:,i);
    if sum(mod(acol,1),'all') ~= 0
        asym = sym(acol);
        [numerators(:,1),denominators(:,1)] = numden(asym);
        commonMultiple = prod(denominators);
        acol = acol*commonMultiple;
        GCD = gcd(acol);
        acol = acol/GCD;
        a(:,i) = acol;
        
    end
end

if [0 1] == size(qi)
    syms q [n 1]; % number of qs ex: L M and T
elseif [n 1] == size(qi)
    q = str2sym(qi);
else
    throw(MException("joshBuckPiTheory:invalidInput","qi is wrong size.")) 
end

syms pi_vecs [n p] % pi's by terms, not a product yet
for i = 1:p
    pi_vecs(:,i) = q.^a(:,i);
end
syms PI [1 p]
for i = 1:p
    PI(i) = prod(pi_vecs(:,i)); % product generates p nondimentional params pi(1-p)
end
end

