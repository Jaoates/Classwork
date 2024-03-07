function [Xdot] = SC_ODE(t,X,I,params,GeoProperties)
% Tfb = FEEDBACK_ODE(X,params); We turned off controls :(
Tfb = [0;0;0];
Tg = gravgrad(X,I);
Ts = SolarPressureTorque(X,GeoProperties);
Ta = AeroTorque(X,GeoProperties);
Tm = MagTorque(t,X,GeoProperties);

T = Tfb+Tg+Ts+Ta+Tm;
Xdot = EOM_ODE(t,X,I,T);


end

