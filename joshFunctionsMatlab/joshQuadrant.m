function [quadrantInt,theta2,quadrantDec] = joshQuadrant(theta1)
% values in rads
% for any theta1 it will return a corresponding positive theta2 as well as
% the quadrant it is in and how many quadrants it is into the unit circle.
% ie, floor of quadrantDec is the quadrant its in.
arguments
    theta1 (1,1) double {mustBeReal,mustBeNonNan}
end
s = sign(theta1);
if s == 0
    theta2 = 0;
    quadrantInt = 0;
    quadrantDec = 0;
    return
end
    theta2 = mod(theta1,2*pi);
    quadrantDec = theta2/(pi/2)+1;
if mod(quadrantDec,1) == 0
    quadrantInt = 0;
else
    quadrantInt = floor(quadrantDec);
end