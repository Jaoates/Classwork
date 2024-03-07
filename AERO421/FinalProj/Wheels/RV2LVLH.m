function [q,eta,epsilon] = RV2LVLH(R,V)
% Takes R and V vectors in ECI and returns te quaternion that roates from ECI to
% LVLH

G.x = [1 0 0]';
G.y = [0 1 0]';
G.z = [0 0 1]';

L.z = -R/norm(R);
h = cross(R,V);
L.y = -h/norm(h);
L.x = cross(L.y,L.z);
joshIsRHONB([L.x,L.y,L.z]);
q = dcm2quat(basis2dcm(L,G));
eta = q(1);
epsilon = q(2:4)';


end


function C21 = basis2dcm(b2,b1)
    C21 = [
    dot(b2.x,b1.x) dot(b2.x,b1.y) dot(b2.x,b1.z)
    dot(b2.y,b1.x) dot(b2.y,b1.y) dot(b2.y,b1.z)
    dot(b2.z,b1.x) dot(b2.z,b1.y) dot(b2.z,b1.z)
    ];
end