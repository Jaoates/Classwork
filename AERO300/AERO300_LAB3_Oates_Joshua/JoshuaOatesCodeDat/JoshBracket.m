function [output] = JoshBracket(f, x0, stepSize, numRoots, searchDomain)

% Joshua Oates - Prelab 3 - algorithim for braketing function
 
% This function will search a domain for roots.
% If no domain is given it will search a certain number of steps up and down
% from an initial guess.
% If it finds the expected number of roots in the initial guesses, it will
% quit early and return all found roots.
% The search will being at the bottom of the domain.
% This does not really represent an ineffecientcy since we are assuming that
% there are multiple roots and thereby the guess is not actually very
% enlightening.

% Note: it will search the entire search domain (inclusive top and bottom) even if the range is not
% divisible by the step size
% Note: it will search nothing outside of the search domain 
% Note: brakets are inclusive on the bottom
% Note: to search an expanding range, this function will have to be used
% recurrsively
% Note: this function will not detect roots which cause the function to
% change sign for less than one stepSize
% Note: this may not detect roots that have a max/min on x axis
% (multipicity even)
arguments
    f
    x0 (1,1) {mustBeNumeric,mustBeReal,mustBeNonNan}
    stepSize (1,1) {mustBeNumeric,mustBeReal,mustBeNonNan}
    numRoots (1,1) {mustBeNumeric,mustBeReal,mustBeNonNan}
    searchDomain (1,2) {mustBeNumeric,mustBeReal,mustBeNonNan} = [0,0]
% get f(x) , x0 , stepSize , numRoots , [searchDomain]
end
 
% if searchDomain DNE %create default searchDomain
%     searchDomain = x0 (plus and minus) stepSize*1000
% end 
if searchDomain(1) == searchDomain(2)
    searchDomain(1) = x0 - stepSize*1000;
    searchDomain(2) = x0 + stepSize*1000;
end
% check inputs are valid
% errID = 'myComponent:inputError';
% msgtext = 'Input does not have the expected format.';
% 
% ME = MException(errID,msgtext)

if (x0 > searchDomain(2)) or (x0 < searchDomain(1))
    error("Initial guess is not in Search Domain")
end




start = searchDomain(1);
lastSign = sign(f(start));
if lastSign == -1
    lastSign = 0;
end

roots = [];
range = abs(searchDomain(1)-searchDomain(2));
I = round((range/stepSize)-mod((range/stepSize),1)); % (Integer for number of steps)

for i = 1:1:I
    thisSign = sign(f(start + stepSize*i));
    if thisSign == -1
        thisSign = 0;
    end
    if length(roots) == numRoots
        break
    end
 
    if thisSign ~= lastSign
        roots = [roots,i];
        lastSign = thisSign;
    end
end

output = [];
for i = 1:1:length(roots)
    r = roots(i);
    output(i,1) = (start + stepSize*(r-1));
    output(i,2) = (start + stepSize*r);
end

thisSign = sign(f(searchDomain(2)));
    if thisSign == -1
        thisSign = 0;
    end
if thisSign ~= lastSign
    i=i+1;
    output(i,1) = (start + stepSize*I);
    output(i,2) = searchDomain(2);
end
