function [theta,M,E] = joshAnomalyCalculator(ecc,anomaly,input)
% M will be Me,Mp or Mh depending on ecc
% E will be Eccentric Anomoly when applicable or F: hyperbolic Ecctric
% anomoly. E will be set to
% values in Rads
arguments
    ecc (1,1) double {mustBeReal,mustBeNonnegative,mustBeNonNan}
    anomaly (1,1) double {mustBeReal,mustBeNonNan}
    input {mustBeMember(input,{'theta','Me'})} = 'theta'
end

if strcmp(input,'theta')
    theta = anomaly;
    if ecc <1 % Me & E
        E = 2*atan(sqrt((1-ecc)/(1+ecc))*tan(theta/2)); % definintion of E, rewriten to solve E
        M = E-ecc*sin(E); % definition of M
    elseif ecc > 1 % Mh & F
        E = log((sqrt(ecc+1)+sqrt(ecc-1)*tan(theta/2))/(sqrt(ecc+1)-sqrt(ecc-1)*tan(theta/2)));
        M = ecc(sinh(F)-F);
    else % ecc == 1 Mp
        E = nan; % This is a rare case and E doesnt have a definition for ecc == 1
        M = .5*tan(theta/2)+(1/6)*tan(theta/2)^3;
    end
elseif strcmp(input,'Me')
    if ecc>=1
        throw(MException("joshAnomalyCalculator:inputNotSupported","entering Me for non eliptical orbits is not supported"))
    end

    M = anomaly;
    if M<pi % initial E
        E0 = M+ecc/2;
    else
        E0 = M-ecc/2;
    end
    f = @(E) (E-ecc*sin(E)-M);
    E = fzero(f,E0);
    theta = 2*atan( tan(E/2)/sqrt((1-ecc)/(1+ecc)) );
end


end