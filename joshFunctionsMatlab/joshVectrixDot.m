function [C] = joshVectrixDot(F1,F2)
% takes the vectrix dot product where the rotation matrix C21 =
% joshVectrixDot(F1,F2) such that r represented in F2 can be found from 
% r2 = C21 * r1. F1 and F2 are 3x3 matricies with columns [x,y,z] where x y
% and z are column vectors representing the basis vectors of F1 and F2

x1 = F1(:,1);
y1 = F1(:,2);
z1 = F1(:,3);

x2 = F2(:,1);
y2 = F2(:,2);
z2 = F2(:,3);

C = [
    dot(x2,x1) dot(x2,y1) dot(x2,z1)
    dot(y2,x1) dot(y2,y1) dot(y2,z1)
    dot(z2,x1) dot(z2,y1) dot(z2,z1)];
end