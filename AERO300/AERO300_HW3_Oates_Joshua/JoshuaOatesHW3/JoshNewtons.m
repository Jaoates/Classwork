function [r, count, xVector, errorVector, errorRatioVector] = JoshNewtons(f, fp, x0, TOL, maxI)
% this function will check if newtons method will converge
% this function takes a diff'able function and its derivative as well as a
% guess and a tolerance.
% this function is not gaurranted to converge and will exit early if it
% does not returning an empty array for r

% this algorithim taken from canvas has been cleaned up for readablity and
% commented for understanding

arguments
    f 
    fp
    x0 (1,1) {mustBeNumeric,mustBeReal}
    TOL (1,1) {mustBeNumeric,mustBeReal} = .0001
    maxI (1,1) {mustBeNumeric,mustBeReal} = 200
end

x1 = x0 - f(x0)/fp(x0); % create x1 for use in algorithim
count = 1;
xVector = [x0; x1];
error = abs(x1 - x0);
errorVector = error;
errorRatioVector = [];

while (error > TOL) && (count <= maxI)
    % set x0 to x1 and recalculate x1 using function and functionPrime
    x0 = x1;
    x1 = x1 - f(x0)/fp(x0); % newtons algorithim
    
    count = count + 1;
    
    % for the sake of the error ratio compute errorSquared
    errorSquared = error^2; % this is the previous error since it hasnt been updated
    error = abs(x1 - x0); % update error 
    errorRatio = error/errorSquared; % update errorRatio
    
    xVector = [xVector; x1]; % append new data to old 
    errorVector = [errorVector; error];
    errorRatioVector = [errorRatioVector; errorRatio];
end

if (count > maxI)
    r = []; % r DNE for these inputs
    warn = "Did not converge in " + num2str(maxI) + " intervals";
    warning(warn)
else
    r = xVector(end);% last guess is also best prediction for the root
end