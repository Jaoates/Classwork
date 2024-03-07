function [r,v,f,f_dot,g_dot] = joshChi2State(X,Cc,Sc)

warning("This hasnt been tested or completed")

f = 1-(X^2/norm(r0))*C(X^2/a);
g = dt-((1/sqrt(mu_e))*X^3*S(X^2/a));
r = f*r0+g*v0;
f_dot = (sqrt(mu_e)/(norm(r)*norm(r0)))*X*((X^2/a)*S(X^2/a)-1);
g_dot = 1-(X^2/norm(r))*C(X^2/a);
v = f_dot*r0 + g_dot*v0;
end

