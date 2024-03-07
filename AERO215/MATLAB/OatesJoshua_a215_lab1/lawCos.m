function [sideC] = lawCos(sideA,sideB,angleAB)
%Uses (sideA,sideB,angleAB) to calculate sideC
%angleAB in degrees
%this is the law of cosines
sideC = sqrt(sideA^2 + sideB^2 - (2*sideA*sideB*cosd(angleAB)));
end


