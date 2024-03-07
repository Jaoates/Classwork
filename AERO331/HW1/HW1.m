%% housecleaning
clear all
close all
clc

addpath('C:\joshFunctionsMatlab\')
%% Problem 2


sympref('FloatingPointOutput',true)
% sympref('FloatingPointOutput',false)

% is the strain tensor
syms E [3 3]
% k is a postive scalar from problem statment
syms k
assume(k>0)

% info on strain state from problem statement
E(3,3) = 0; % Ezz
E(1,3) = 0; % Exz
E(2,3) = 0; % Eyz

% strain gauges read this in problem statement
E1 = -k;
E2 = k;
E3 = -k/2;

% symetric condition
E(2,1) = E(1,2);
E(3,1) = E(1,3);
E(3,2) = E(2,3);


% direction vectors of strain guages
n1 = [cos(pi/6);sin(pi/6);0];
n2 = [1;0;0];
n3 = [cos(pi/6);-sin(pi/6);0];

% systems of eqn
eqn1 = E1 == n1'*E*n1;
eqn2 = E2 == n2'*E*n2;
eqn3 = E3 == n3'*E*n3;

% solve for Eyy and Exz and Exx
sol = solve([eqn1,eqn2,eqn3],[E1_2,E2_2,E1_1]);
disp(sol)

% put them in the Matrix
E(1,1) = sol.E1_1;
E(1,2) = sol.E1_2;
E(2,2) = sol.E2_2;

% symetric condition
E(2,1) = E(1,2);
E(3,1) = E(1,3);
E(3,2) = E(2,3);

% disp answer
simplify(E)

% rename E1-3 to be principle stains, instead of input strains
clear E1 E2 E3 eqn1 eqn2 eqn3

% evecs, Ep is E'
[evecs,Ep] = eig(E);

[V,D] = joshBasisFix(evecs,Ep)

% syms lam
% syms e1 e2 e3 [3 1]
% 
% eqn = det(E-lam*eye(3)) == 0
% sol = simplify(solve(eqn))
% 
% lam = sol(1);
% eqn = (E-lam*eye(3))*e1 == 0;
% e1 = simplify(solve(eqn,e1))
% 
% lam = sol(2);
% eqn = (E-lam*eye(3))*e2 == 0;
% e2 = simplify(solve(eqn,e2))
% 
% lam = sol(3);
% eqn = (E-lam*eye(3))*e3 == 0;
% e3 = simplify(solve(eqn,e3))

Ep = simplify(Ep);
% e1 is i' etc
e1 = evecs(:,1);
e2 = evecs(:,2);
e3 = evecs(:,3);

% Ep

Ep = [[-Ep(2,2),0,0];...
       [0,-Ep(1,1),0];...
       [0,0,-Ep(3,3)]];
Ep = [Ep(1,1),Ep(2,2),Ep(3,3)];
e1 = -evecs(:,2); % x becomes -y
e2 = -evecs(:,1); % y becomse -x
e3 = -evecs(:,3); % z flips
% clear evecs
Epsave = [
[Ep(1),0,0]    
[0,Ep(2),0]
[0,0,Ep(3)]
];

e1 = simplify(e1/norm(e1))
e2 = simplify(e2/norm(e2))
e3 = simplify(e3/norm(e3))

evecs(:,1) = e1
evecs(:,2) = e2
evecs(:,3) = e3

% plug in k for numeric answer
Epk = subs(Ep,k,.002);

gamMax = Ep(1)-Ep(3)
clear n1 n2 n3
n1 = (sqrt(2)/2)*(e1+e3)
n2 = (sqrt(2)/2)*(e1-e3)

%% Part 2
