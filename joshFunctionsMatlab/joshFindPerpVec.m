function [vb1,vb2,vb3] = joshFindPerpVec(v)
% this function will find 2 vectors which are perpendicular to eachother
% and to the input vector v. They are not any particular perpendicular
% vectors. vp3 is the unit vector of the original vector v. vp1,vp2 and vp3
% will for an RHONB.

% outBasis z
vb3 = v/norm(v);

% arbitrary basis
b = eye(3);
c = [cross(vb3,b(:,1)),cross(vb3,b(:,2)),cross(vb3,b(:,3))];
[~,I] = max(vecnorm(c)); % take the largest cross and use this as vb1
c = c(:,I);

% RHONB rules
vb1 = c/norm(c); % normalize this vector 
vb2 = -cross(vb1,vb3);
end