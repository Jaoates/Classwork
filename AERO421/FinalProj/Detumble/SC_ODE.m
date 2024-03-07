function [Xdot] = SC_ODE(t,X,I,k)
T = FEEDBACK_ODE(X,k);
Xdot = EOM_ODE(t,X,I,T);
end

