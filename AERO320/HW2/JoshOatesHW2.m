
%% Problem 1 abcd
clear all;
clc

addpath('C:\joshFunctionsMatlab\')

C21 = [[0 0 1];[1 0 0];[0 1 0]]; % Given
myString1 = string(joshIsOnes(C21*C21' == eye(3) & C21'*C21 == eye(3) & det(C21) == 1)); % 1a

[a_v,phi] = joshRotM2PrincAxe(C21) % a_v and phi

phi_sym = sym(phi); % cast symbolic for readablility
a_v_sym = sym(a_v);

C21_star = joshPrincAxe2RotM(a_v,-phi);
C21_star = round(C21_star,15) % round off any errors near e-mach

myString2 = string(isequal(C21',C21_star)); % C21' == C21_star

C21_pound = joshPrincAxe2RotM(-a_v,-phi);
C21_pound = round(C21_pound,15)

myString3 = string(isequal(C21,C21_pound)); % C21 == C21_pound

disp("My workings for Problem 1 have the following results:")
disp("C21 is a rotation matrix: "+myString1)

disp("axis vector a is: ")
disp("   " + string(a_v_sym))
disp("rotation angle phi is: "+string(phi_sym))

disp("C21 transposed is = to C21*: "+myString2)
disp("C21 is = to C21#: " + myString3)


%% Problem 2
clear all;

Cx = @(theta)...
    [[1 0 0];...
    [0 cos(theta) sin(theta)];...
    [0 -sin(theta) cos(theta)]];

Cy = @(theta)...
    [[cos(theta) 0 -sin(theta)];...
    [ 0 1 0];...
    [sin(theta) 0 cos(theta)]];

Cz = @(theta)...
    [[cos(theta) sin(theta) 0];...
    [-sin(theta) cos(theta) 0];...
    [0 0 1]];

syms t  [3 1] % x y and z thetas
assume(t,'real');

Cy1 = Cy(t(1)) % generate the individual rotation matrices
Cz1 = Cz(t(2))
Cx1 = Cx(t(3))

C21 = Cx1*Cz1*Cy1;

C21s1 = subs(C21,t(1),pi/2)
C21s2 = subs(C21,t(2),pi/2)
C21s3 = subs(C21,t(3),pi/2)

disp("My workings for Problem 2 have the following results:")
disp("C21 for a 2-3-1 rotation is given by: ")
disp("   "+string(C21))
disp("By analysis of the matrices resulting from substituting pi/2 for theta-x, theta-y, theta-z respectively, only theta-z = pi/2 is not solveable. The singularity of 2-3-1 is theta-z.")


%% Problem 3 part 1 b
clear all;

syms a [3 1]
assume(a,'real') % is vector of real numbers
assumeAlso(sqrt(sum(a.^2))==1) % is unit vector
% assumptions(a)
ax = joshCross(a)
LHS = ax*ax*ax
RHS = -ax
myString = string(joshIsOnes(isAlways(RHS == LHS))); % returns a matrix of logical and checks each individual LHS value is always the corresponding value on RHS

disp("My workings for Problem 3 part 1 b have the following results:")
disp("axaxax is = to -ax: "+myString)

%% Problem 3 part 2
clear all;

Cx = @(theta)...
    [[1 0 0];...
    [0 cosd(theta) sind(theta)];...
    [0 -sind(theta) cosd(theta)]];

syms A B 

LHS = Cx(A+B)
RHS = Cx(A)*Cx(B)
myString = string(joshIsOnes(isAlways(RHS == LHS)));
disp("My workings for Problem 3 part 2 have the following results:")
disp("Cx(A+B) is = to Cx(A)+Cx(B): "+myString)


