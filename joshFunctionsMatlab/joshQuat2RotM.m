function [C1,C2] = joshQuat2RotM(eta,epsilon)
C2 = (2*eta^2-1)*eye(3)+2*epsilon*epsilon'-2*eta*joshCross(epsilon);
end