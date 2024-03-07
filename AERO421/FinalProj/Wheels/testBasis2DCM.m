

addpath("C:\joshFunctionsMatlab\")

v = [1 1 0]';
[x1,y1,z1]=joshFindPerpVec(v);
b1.x = x1;
b1.y = y1;
b1.z = z1;

G.x = [1 0 0]';
G.y = [0 1 0]';
G.z = [0 0 1]';

basis2dcm(b1,G)

disp("end")
function C21 = basis2dcm(b2,b1)
    C21 = [
    dot(b2.x,b1.x) dot(b2.x,b1.y) dot(b2.x,b1.z)
    dot(b2.y,b1.x) dot(b2.y,b1.y) dot(b2.y,b1.z)
    dot(b2.z,b1.x) dot(b2.z,b1.y) dot(b2.z,b1.z)
    ];
end