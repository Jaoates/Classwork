function [cirPts] = joshRandCircle(n)
% returns n number of points on the unit circle as column vectors.
arguments
    n = 1
end

theta = rand([1,n])*2*pi;
cirPts = [cos(theta);sin(theta)];
end