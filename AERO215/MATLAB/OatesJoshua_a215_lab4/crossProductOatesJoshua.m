function [Out] = crossProductOatesJoshua(H1,H2)
%   Takes two 3x1 vectors and computes cross product
Out(1) = H1(2)*H2(3) - H1(3)*H2(2);
Out(2) = H1(3)*H2(1) - H1(1)*H2(3);
Out(3) = H1(1)*H2(2) - H1(2)*H2(1);

end

