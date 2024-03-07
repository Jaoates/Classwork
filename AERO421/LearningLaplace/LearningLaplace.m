clear all
clc

syms x(t) t
dx = diff(x)
ddx = diff(dx)

b = [1 -2]

f = ddx+b(1)*dx+b(2)*x == 0

simplify(ilaplace(laplace(f)))


