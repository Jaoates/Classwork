function [eps, eta] = rot2quat(C)
    eta = .5*(1+trace(C))^.5;
    
    epsx = (C(2,3)-C(3,2))/(4*eta);
    epsy = (C(3,1)-C(1,3))/(4*eta);
    epsz = (C(1,2)-C(2,1))/(4*eta);
    eps = [epsx; epsy; epsz];
end