%% 0 
close all;
clear all;
clc
addpath("C:\joshFunctionsMatlab\")  
%% Problem 1
[Cx,Cy,Cz] = joshAxisRotation();

%a


r_1= [6738;3391;1953];
v_1=[-3.5;4.39;4.44];
zl_1 = r_1/norm(r_1);
yl_1 = cross(r_1,v_1)/norm(cross(r_1,v_1));
xl_1 = cross(yl_1,zl_1);
xe_1 = [1;0;0];
ye_1 = [0;1;0];
ze_1 = [0;0;1];

C21 = [[dot(xl_1,xe_1),dot(xl_1,ye_1),dot(xl_1,ze_1)];[dot(yl_1,xe_1),dot(yl_1,ye_1),dot(yl_1,ze_1)];[dot(zl_1,xe_1),dot(zl_1,ye_1),dot(zl_1,ze_1)]];

%% bcd
[a,phi] = joshRotM2PrincAxe(C21);

thetax = atan2(C21(2,3),C21(3,3));
thetay = asin(-C21(1,3));
thetaz = atan2(C21(1,2),C21(1,1));

C21verify = Cx(thetax)*Cy(thetay)*Cz(thetaz); % for my use

eta =@(C) .5*sqrt(1+trace(C));
epsilon =@(C,eta)...
    [(C(2,3)-C(3,2))/(4*eta);...
    (C(3,1)-C(1,3))/(4*eta);...
    (C(1,2)-C(2,1))/(4*eta)];
eta = eta(C21);
epsilon = epsilon(C21,eta);
disp("------------P1------------")
disp("My workings for this problem have the following results:")
disp("Rotation matrix from ECI to the spacecraft: ")
disp(round(C21,2))
disp("Principle axis (a) for this transformation: ")
disp(round(a,2))
disp("Principle angle of this transformation: "+ string(phi))
disp("Euler angles:")
disp("      thetax: " + string(round(rad2deg(thetax),2))+" degrees")
disp("      thetay: " + string(round(rad2deg(thetay),2))+" degrees")
disp("      thetaz: " + string(round(rad2deg(thetaz),2))+" degrees")
disp("quaternion:")
disp("      eta: "+string(round(eta,2)))
disp("      epsilon: ")
disp("      "+round(epsilon,2))


%% Problem 2
clear all;
[Cx,Cy,Cz] = joshAxisRotation();

syms t [3 1]

C21oft = Cz(t(3))*Cx(t(2))*Cz(t(1));

clear t t1 t2 t3

syms C [3 3]

t2 = acos(C(3,3)); %theta
t1 = asin(C(3,1)/sin(t2)); %phi
t3 = acos(C(1,1)/cos(t1)); %psi

disp("------------P2------------")
disp("My workings for this problem have the following results:")
disp("The formula for finding C21 in terms of phi = t1, theta = t2, and psi = t3, is:")
disp(C21oft)
disp("the formulas for finding the Euler angles are: ")
disp("t2 = "+string(t2))
disp("t1 = "+string(t1))
disp("t3 = "+string(t3))
disp("The singularities for this scheme are at t2 = 0 and t2 = pi. Physically we can consider a plane which rolls pitches and rolls again. If the plane ends its rotation pointing the oposite direction from where it began (a pitch of pi) or continues to point in the original direction (a pitch of 0) then the last degree of freedom to define is the roll of the plane. This can be defined by either the first or second roll. How much the first vs second roll contributed cannot be determined from soley the initial and final orientations of the plane.")
%% Problem 3
clear all;
[Cx,Cy,Cz] = joshAxisRotation();

n= sqrt(2)/2
C21=[[n 0 -n];[0 1 0];[n 0 n]]
C32=[[0 0 -1];[-1 0 0];[0 1 0]]
clear n

%% a
isRotM = @(M) (round(M*M',14) == eye(3) & round(M'*M,14) == eye(3) & round(det(M),14) == 1)

isRot21 = joshIsOnes(isRotM(C21));
isRot32 = joshIsOnes(isRotM(C32));
%% b
eta =@(C) .5*sqrt(1+trace(C));
epsilon =@(C)...
    [(C(2,3)-C(3,2))/(4*eta(C));...
    (C(3,1)-C(1,3))/(4*eta(C));...
    (C(1,2)-C(2,1))/(4*eta(C))];

C=@(eta,epsilon) (2*eta^2-1)*eye(3)+2*epsilon*epsilon'-2*eta*joshCross(epsilon)

eta21 = eta(C21);
epsilon21 = epsilon(C21);
eta32 = eta(C32);
epsilon32 = epsilon(C32);
C21verify = C(eta21,epsilon21);% for my use 
C32verify = C(eta32,epsilon32);

%% c
eta31=eta21*eta32-epsilon32'*epsilon21;
epsilon31 = eta32*epsilon21+eta21*epsilon32+joshCross(epsilon21)*epsilon32;

%% de
C31 =C32*C21;
eta31verify = eta(C31);
epsilon31verify = epsilon(C31);
isSame1 = round(eta31verify,14) == round(eta31,14);
isSame2 = round(epsilon31verify,14) == round(epsilon31,14)
isSame2 = joshIsOnes(isSame2)

%% f
eta21inv = eta(inv(C21));
epsilon21inv = epsilon(inv(C21));

isSame3 = round(eta21inv,14) == round(eta21,14);
isSame4 = round(epsilon21inv,14) == round(-epsilon21,14)
isSame4 = joshIsOnes(isSame4)

disp("------------P3------------")
disp("My workings for this problem have the following results:")
disp("C21 is a rotation matrix: "+string(isRot21))
disp("C32 is a rotation matrix: "+string(isRot32))
disp("quaternion21:")
disp("      eta: "+string(eta21))
disp("      epsilon: ")
disp("      "+epsilon21)
disp("quaternion32:")
disp("      eta: "+string(eta32))
disp("      epsilon: ")
disp("      "+epsilon32)
disp("quaternion31:")
disp("      eta: "+string(eta31))
disp("      epsilon: ")
disp("      "+epsilon31)
disp("C31:")
disp(string(C31))
disp("quaternion31 from C31:")
disp("      eta: "+string(eta31verify))
disp("      epsilon: ")
disp("      "+epsilon31verify)
disp("eta31 has been found both ways: "+string(isSame1))
disp("epsilon31 has been found both ways: "+string(isSame2))
disp("eta21 == eta21^-1: "+string(isSame3))
disp("-epsilon21 == epsilon21^-1: "+string(isSame4))
disp("The quaternions are related using:")
disp("-epsilon21 == epsilon21^-1, and eta21 == eta21^-1")

