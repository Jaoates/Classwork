function  [e_vec, e_norm, Te, T] = geo_monte_sim(r0, v0, T0lla, jd, err, n)

for i = 1:n
[T(:,i),Te(:,i)] = geo_pointing(r0,v0,T0lla,jd,err);
end
e_vec = T-Te;
e_norm = vecnorm(e_vec);
end