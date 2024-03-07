function [varargout] = joshInterp(x1,y1,x2,y2,x)
% [y] = joshInterp(x1,y1,x2,y2,x)
% [S] = joshInterp(x1,y1,x2,y2)
% [y,S] = joshInterp(x1,y1,x2,y2,x)


% you can do this by hand or with a poly fit, but this is easy
% takes 2 points and linearly interpolates or extrapolates to give you a
% value y corresponding to x

% overloaded here so that if x is not provided, y will be a structure with 
arguments
x1 (1,1) {mustBeReal, mustBeNumeric}
y1 (1,1) {mustBeReal, mustBeNumeric}
x2 (1,1) {mustBeReal, mustBeNumeric}
y2 (1,1) {mustBeReal, mustBeNumeric}
x {mustBeReal, mustBeNumeric,mustBeVector} = nan
end

if x2 <= x1
    throw(MException('joshInterp:invalidInput','x2 must be greater than x1'))
end

if nargout > 2
    throw(MException('joshInterp:invalidOutput','must return either 1 or 2 output arguments'))
end

% logic
m = (y2-y1)/(x2-x1);
b = y1-m*x1;

y = m*x+b; % might end up y is nan if x is nan, y shouldn't return if y is nan

% prep for output
myStruct.m = m;
myStruct.yInt = b;
myStruct.xInt = -b/m;

% output
varargout = cell(nargout,1);

if isnan(x)
    varargout{1} = myStruct;
elseif nargout() == 2
    varargout{1} = y;
    varargout{2} = myStruct;
else
    varargout{1} = y;
end

end