function [eta,epsilon] = joshRotM2Quat(C)

% functions works if eta ~= 0, this will be patched later
if ~joshIsRotM(C,1e-4)
    throw(MException("joshRotM2Quat:invalidInput","C must be a rotational matrix"))
end

eta =.5*sqrt(1+trace(C));
epsilon =...
    [(C(2,3)-C(3,2))/(4*eta);...
    (C(3,1)-C(1,3))/(4*eta);...
    (C(1,2)-C(2,1))/(4*eta)];
end

