function [phiA,phiD,Tsyn,Twait,Tlw] = joshPhaseAngle(r1,r2,T,mu)
% returns phase angle and synodic period

% r1 and r2 as magnitudes
% T is transfer time, same from joshHomann

% from D to A (depart and arrive)



arguments
    r1
    r2
    T
    mu
end
warning("joshPhaseAngle: argument validation is incomplete, function is untested")
warning("joshPhaseAngle: assume only HT with circular orbit")

n1 = sqrt(mu/r1^3);
n2 = sqrt(mu/r2^3);

T1 = (2*pi/sqrt(mu))*r1^1.5;
T2 = (2*pi/sqrt(mu))*r2^1.5;
Tsyn = (T1*T2)/abs(T1-T2); % synotic period

% phiA = thetaA - thetaD - n1*T;
phiA = pi - n1*T;
phiD = pi - n2*T;

% N is a positive integer, It should be sufficently large so that Twait is positive, I think we can solve this with a mod...
if n1>n2
    % Twait= (-2*phiA - (2*pi)*N)/(n2-n1);
    Twait= (-2*phiA - (2*pi))/(n2-n1); % how long until I can transfer back?
else
    Twait= (-2*phiA + (2*pi))/(n2-n1);
end
Twait = mod(Twait,2*pi);


% similar to N above, gets you into the window/s you want
% Tepoch is like J2000?
Tlw = (phiD-phiA/(n2-n1))+Tepoch+Tsyn+numTsyn; % launch window

end