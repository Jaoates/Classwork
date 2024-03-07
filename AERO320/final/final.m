clear all
close all
clc
syms x r mg l
diff(sin(x))

r=[-r,0,0]
f=[0,-mg,0]

cross(r,f)

eq = mg*l*(1-cos(x))*cos(x)
deq = diff(eq,x)
diff(cos(x))
diff(cos(x)^2) 
isAlways(deq == mg*l*(2*cos(x)*sin(x)-sin(x)))
deq = mg*l*(2*cos(x)*sin(x)-sin(x))
assume(-2*pi<x<2*pi)
assumeAlso(x,'real')
isAlways(deq<=0)

diff(cos(x)*cos(x))