a2 =   108600000;
a1 =   149600000;
mu =   1.3270e+11;
n1 =    0.9856 % deg/day
n2 =    1.6020
at = (a2+a1)/2


Ttrans = (2*pi/sqrt(mu))*at^1.5;
Ttrans = Ttrans/(60*60*24)
phiD = 180-(n2*Ttrans)
phiA = 180-(n1*Ttrans)

% 1960 jan1.5
tepoch  = 2436934.5;%days
theta1 = 174.2943; % location at epoch
theta2 = 100.1582;

Tsyn = 583.92;
n = ceil((tlaunchdate - tepoch)*365.25/Tsyn) % this is how many periods have passed from epoch to today, 

% answer: 2443734: aug 14th, 1978
