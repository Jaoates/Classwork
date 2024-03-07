

%% 14.1

syms Ix Iy Iz nu k ks

A = [-k/Ix -(Iz-Iy)*nu/Ix
    -(Ix-Iz)*nu/Iy -k/Iy ]

assume(Ix>Iz)
assumeAlso(Iz>Iy)
assumeAlso(Iy>0)
assumeAlso(nu>0)
assumeAlso(k>0)

pole = simplify(eig(A))
% pole = subs(pole,k,sqrt(nu^2*(Iz-Iy)*(Ix-Iz)))
% isAlways(pole<=0)

% copied output for expediency
-(Ix*k + Iy*k + (- 4*Ix^2*Iy^2*nu^2 + 4*Ix^2*Iy*Iz*nu^2 + Ix^2*k^2 + 4*Ix*Iy^2*Iz*nu^2 - 4*Ix*Iy*Iz^2*nu^2 - 2*Ix*Iy*k^2 + Iy^2*k^2)^(1/2))/(2*Ix*Iy)
-(Ix*k + Iy*k - (- 4*Ix^2*Iy^2*nu^2 + 4*Ix^2*Iy*Iz*nu^2 + Ix^2*k^2 + 4*Ix*Iy^2*Iz*nu^2 - 4*Ix*Iy*Iz^2*nu^2 - 2*Ix*Iy*k^2 + Iy^2*k^2)^(1/2))/(2*Ix*Iy)

a = - 4*Ix^2*Iy^2*nu^2 + 4*Ix^2*Iy*Iz*nu^2 + Ix^2*k^2 + 4*Ix*Iy^2*Iz*nu^2 - 4*Ix*Iy*Iz^2*nu^2 - 2*Ix*Iy*k^2 + Iy^2*k^2
solve(a == ((Ix+Iy)*k)^2,k)

% for ks

%% 15.2


