function [int1,int2] = joshSphereIntersect(u,c,r)
% finds the intersection(s) of a sphere and the line which is defined by
% the vector u. int1 and 2 are the location of the intersetion relative to
% the vectors origin.

% c is the center of the sphere relative to the orgin of the vector
% r is the radius of the sphere (r = 1 if not specified)

% o is a point through which the line passes

%https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

arguments
    u
    c
    r = 1
end

o = [0;0;0]; % set it so the vector is always at the orign
u = u/norm(u);
del = dot(u,(o-c))^2-(norm(o-c)^2-r^2);
d1 = -(dot(u,(o-c)))-sqrt(del);
d2 = -(dot(u,(o-c)))+sqrt(del);

int1 = u*d1;
int2 = u*d2;


end