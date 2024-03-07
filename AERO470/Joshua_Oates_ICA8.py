import numpy as np
import matplotlib.pyplot as plt
from vpython import *
# from scipy import linalg as spl
# from scipy import integrate as int

sun = sphere(pos=vec(0,0,0), radius=0.4, color=color.yellow)
earth = sphere(pos=vec(0,0,0), radius=0.1, color=color.blue)
moon = sphere(pos=vec(0,0,0), radius=0.05, color=color.white)

def circularMotion(t,r,v):
    return vec(cos(t*v),sin(t*v),0)*r

t = 0
dt = .01
re = 3
rm = .4
ve = 1
vm = 4

tend = 5
while t<tend:
    rate(30)
    earth.pos = circularMotion(t,re,ve)
    moon.pos = earth.pos + circularMotion(t,rm,vm)
    t = t+dt
    
t = 0
rc = 4
vc = .3

plane = box(pos=vec(0,0,0),length = 1,height = 1,width = .1)

while 1:
    rate(30)
    earth.pos = circularMotion(t,re,ve)
    moon.pos = earth.pos + circularMotion(t,rm,vm)

    scene.camera.pos = vec(0,6,sin(t*3)*3)
    scene.camera.up = vec(0,0,1)
    scene.camera.axis = -scene.camera.pos

    t = t+dt

input("End of line [ENTER]")

