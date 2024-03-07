function [r, count, xVector, errorVector, errorRatioVector] = JoshSecant(f,x0,x1, TOL, maxI)
% this function will check if newtons method will converge
% this function takes a function and two initial guesses as well as a tolerance.
% this function is not gaurranted to converge and will exit early if it
% does not returning an empty array for r

% this algorithim taken from canvas has been cleaned up for readablity and
% commented for understanding

arguments
    f 
    x0 (1,1) {mustBeNumeric,mustBeReal}
    x1 (1,1) {mustBeNumeric,mustBeReal}
    TOL (1,1) {mustBeNumeric,mustBeReal} = .0001
    maxI (1,1) {mustBeNumeric,mustBeReal} = 200
end

if x0 == 0 || x1 == 0
    warning("input guesses should be non zero")
end
if x0 == x1
    warning("input guesses should not be the same value")
end

x2 = x0 - f(x0)*((x1 - x0)/(f(x1)-f(x0))); % create x1 for use in algorithim
count = 1;
xVector = [x0; x1;x2];
error = abs(x1 - x0);
errorVector = error;
errorRatioVector = [];


while (error > TOL) && (count <= maxI)
    % set x0 to x1 and recalculate x1 using function and functionPrime
    x0 = x1;
    x1 = x2;
    x2 = x0 - f(x0)*((x1 - x0)/(f(x1)-f(x0))); % secant algorithim
    
    
    count = count + 1;
    
    % for the sake of the error ratio compute errorSquared
    errorSquared = error^2; % this is the previous error since it hasnt been updated
    error = abs(x2 - x1); % update error 
    errorRatio = error/errorSquared; % update errorRatio
    
    xVector = [xVector; x2]; % append new data to old 
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