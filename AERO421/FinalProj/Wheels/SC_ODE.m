function [Xdot] = SC_ODE(t,X,I,params,GeoProperties)
Mc = FEEDBACK_ODE(X,params);
 
[Mw,Ohmdot] = wheelEOM(X,Mc,GeoProperties);

Tfb = Mw;

Tg = gravgrad(X,I);
Ts = SolarPressureTorque(X,GeoProperties);
Ta = AeroTorque(X,GeoProperties);
Tm = MagTorque(t,X,GeoProperties);

T = Tfb+Tg+Ts+Ta+Tm;

% T = Tfb;

Xdot = EOM_ODE(t,X,I,T);
Xdot = [Xdot;Ohmdot];

end

