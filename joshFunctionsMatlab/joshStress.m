function [varargout] = joshStress(sigma)
sigmaH = (1/3)*trace(sigma);
sigmaDev = sigma - sigmah*eye(3);
sigmaE = ((3/2)*sum(sum(sigmaDev.^2)))^.5;

joshBasisFix
tauMax = 

outStruct.sigmaH = sigmaH;
outStruct.sigmaDev = sigmaDev;
outStruct.sigmaE = sigmaE;
outStruct.tauMax = tauMax;
end