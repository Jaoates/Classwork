function [out] = joshAdvBeam(Ai,yi_prime,zi_prime,Iyoiyoi,Izoizoi,Iyoizoi,Ei_E1,alphai,Ei)
% this function takes several 



% arguments
%     Ai {mustBeReal}
%     yi_prime {mustBeReal}
%     zi_prime {mustBeReal}
%     Iyoiyoi {mustBeReal}
%     Izoizoi {mustBeReal}
%     Iyoizoi {mustBeReal}
%     Ei_E1 {mustBeReal}
%     alphai {mustBeReal} = nan
%     Ei {mustBeReal} = nan
% end

arguments
    Ai 
    yi_prime 
    zi_prime 
    Iyoiyoi 
    Izoizoi
    Iyoizoi
    Ei_E1
    alphai = nan
    Ei  = nan
end

A = Ai;
yp = yi_prime;
zp = zi_prime;
Iz0 = Izoizoi;
Iy0 = Iyoiyoi;
Iyz0 = Iyoizoi;

% n = length(A);
% if length(yp) ~= n  | length(zp) ~= n  | length(Iz0) ~= n  | length(Iy0) ~= n  | length(Iyz0) ~= n | length(Ei_E1) ~= n 
%     throw(MException('joshAdvBeam:invalidInput','At least one of the input vectors is not the correct length'))
% end


% Ai*(Ei/E1)
AE_E1 = Ei_E1.*Ai;
% A*
As = sum(AE_E1);

% A*(E/E1)*y'
AE_E1yp = AE_E1.*yp;
% y'*
yps = sum(AE_E1yp)/As;


% Ai*(Ei/E1)*zi'
AE_E1zp = AE_E1.*zp;
% z'*
zps = sum(AE_E1zp)/As;

% yy
% (Ei/E1)*(Iyoiyoi+Ai'*zi'^2)
var1 = (Ei_E1.*(Iy0+A.*zp.^2));
% I*y'y'
Iyps = sum(var1);
% I*yy = I*y'y' - A*(z'*)^2
Iys = Iyps - As.*zps.^2;

% zz
var2 = (Ei_E1.*(Iz0+A.*yp.^2));
Izps = sum(var2);
Izs = Izps - As.*yps.^2;

% yz
var3 = (Ei_E1.*(Iyz0+A.*zps.*yps));
Iyzps = sum(var3);
Iyzs = Iyzps - As.*zps.*yps;

% y and z
y = yp-yps;
z = zp-zps;

out.y = y;
out.z = z;
out.As = As;
out.yps = yps;
out.zps = zps;

% out.Iyyps = Iyps;
out.Iyys = Iys;

% out.Izzps = Izps;
out.Izzs = Izs;

% out.Iyzps = Iyzps;
out.Iyzs = Iyzs;

if (~isnan(alphai)) & (~isnan(Ei))

    if length(alphai) ~= n | length(Ei) ~= n 
        throw(MException('joshAdvBeam:invalidInput','Either alphai or Ei is the wrong length'))
    end
    E = Ei;

    EalphaA = E.*alphai.*A;
    EalphaAy = E.*alphai.*A.*y;
    EalphaAz = E.*alphai.*A.*z;
    
    PT_DT = sum(EalphaA);
    Mz_DT = sum(EalphaAy);
    My_DT = sum(EalphaAz);

    out.PT_DT = PT_DT;
    out.MzT_DT = Mz_DT;
    out.MyT_DT = My_DT;
end


end