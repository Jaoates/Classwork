function [sideC] = joshLawCos(sideA,sideB,angleAB)
%Uses (sideA,sideB,angleAB) to calculate sideC
%angleAB in rad
%this is the law of cosines
sideC = sqrt(sideA^2 + sideB^2 - (2*sideA*sideB*cos(angleAB)));
end


