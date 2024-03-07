from random import *
from math import *
import matplotlib.pyplot as plt
import numpy as np
#region Q1

class Point():
    def __init__(self,x,y):
        self.x = x
        self.y = y
        self.distance = sqrt(x**2+y**2)
    def __str__(self):
        # return(f"{self.x} , {self.y}")
        return(f"x:{self.x:1.4f} | y:{self.y:1.4f}".format())
    def __repr__(self):
        return(f"Point with {self.__str__()}")
    
def genRandPoint(xdomain=[0,1],ydomain=[0,1]):
    return Point(uniform(xdomain[0],xdomain[1]),uniform(ydomain[0],ydomain[1]))

tol = .01

npts = [0]
inOutRatio = [3]
accuracy = 10
while accuracy > tol:
    npts.append(npts[-1]+10)
    pts = [genRandPoint() for x in range(npts[-1])]
    ptsin = [pt for pt in pts if pt.distance<1]
    # area of a quarter unit circle * 4 is an estimation of pi
    inOutRatio.append(len(ptsin)*4/len(pts))
    # accuracy = abs((sum(inOutRatio)/len(inOutRatio)-inOutRatio[-1]))
    accuracy = abs((pi)-inOutRatio[-1])


print(f"The number of points required to reach an estimation of pi: {inOutRatio[-1]:1.5f} within {tol*100}% is {npts[-1]}")

plt.style.use('_mpl-gallery')

# plot convergence
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.07,left=0.07)
ax.plot(npts, inOutRatio, linewidth=2.0,label='estimate')
ax.hlines(pi,0,npts[-1],label='pi',color = 'red')
fig.legend()

fig.show()

# plot final points graph
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.07,left=0.07)
ax.scatter([pt.x for pt in pts],[pt.y for pt in pts],color = "red")
ax.scatter([pt.x for pt in ptsin],[pt.y for pt in ptsin],color = "blue")
ax.set_aspect("equal","box")

fig.show()
#endregion

#region Q2

def myFunc(x):
    return x**3*sin(x)

domain = np.arange(-2,3,.0001)
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.07,left=0.07)
ax.plot(domain, [myFunc(x) for x in domain], linewidth=2.0,label='function')

area = np.trapz([myFunc(x) for x in domain],domain)
xdomain = [-2,3]
ydomain = [0,10]
boxarea = 50
realRatio = area/boxarea
print(realRatio)
tol = .01

npts = [0]
inOutRatio = [0]
accuracy = 10
while accuracy>tol:
    npts.append(npts[-1]+10)
    pts = [genRandPoint(xdomain,ydomain) for x in range(npts[-1])]
    ptsin = [pt for pt in pts if pt.y<myFunc(pt.x)]
    # area of a quarter unit circle * 4 is an estimation of pi
    inOutRatio.append(len(ptsin)/len(pts))
    # accuracy = abs((sum(inOutRatio)/len(inOutRatio)-inOutRatio[-1]))
    accuracy = abs((realRatio)-inOutRatio[-1])
    print(inOutRatio[-1])

inOutRatio =[x * boxarea for x in inOutRatio]
print(f"The number of points required to reach an estimation of the integral: {inOutRatio[-1]:1.5f} within {tol*100}% is {npts[-1]}")

plt.style.use('_mpl-gallery')

# plot final points graph
# fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.07,left=0.07)
ax.scatter([pt.x for pt in pts],[pt.y for pt in pts],color = "red")
ax.scatter([pt.x for pt in ptsin],[pt.y for pt in ptsin],color = "blue")
ax.set_aspect("equal","box")
fig.show()


# plot convergence
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.07,left=0.07)
ax.plot(npts, inOutRatio, linewidth=2.0,label='estimate')
ax.hlines(area,0,npts[-1],label='real area',color = 'red')
fig.legend()
fig.show()
#endregion

#region Q3
R = 287.053
# rho0 = 1.225
# T0 = 15
T = 25 + 273.15
Ts = .2
P = 104847
Ps = 52

def rho(P,T): return P/(R*T)

n = 1000
rhom=[]
Tm = []
Pm = []
for i in range(n):
    Tm.append(np.random.normal(T,Ts))
    Pm.append(np.random.normal(P,Ps))
    rhom.append(rho(Pm[-1],Tm[-1]))

rhobar = np.mean(rhom)
sig = sqrt(sum((rhom - rhobar)**2)/n)

fig,axs = plt.subplots(2)
fig.subplots_adjust(bottom=0.07,left=0.07)
axs[0].scatter([r-rho(P,T) for r in rhom], [t-T for t in Tm], linewidth=2.0)
axs[0].set_xlabel("del-T")
axs[0].set_ylabel("del-rho")
axs[1].scatter([r-rho(P,T) for r in rhom], [p-P for p in Pm], linewidth=2.0)
axs[1].set_xlabel("del-P")
axs[1].set_ylabel("del-rho")
fig.legend()
fig.show()




print(f"Monte Carlo estimates that the mean of the density is {rhobar:1.4f} with a standard deviation of {sig:1.6f}")
# print(rhobar)
# print(sig)

# partial dervatives of  rho:
drho_dP =1/(R*T)
drho_dT = -P/(R*T**2)
sig = np.linalg.norm([drho_dP*Ps,drho_dT*Ts])

print(f"The analytical method finds that the mean of the density is {rhobar:1.4f} with a standard deviation of {sig:1.6f}")
#endregion

#region Q4


def earnings(p,s,vc,fc): return p*s-vc-fc
def getp(): return triangular(50,70,55)
def gets(): return triangular(2000,3000,2440)
def getvc(): return triangular(50000,65000,55200)
def getfc(): return triangular(40000,125000,65000)

n = 10000
p=[]
s=[]
vc=[]
fc=[]
earned=[]
for i in range(n):
    p.append(getp())
    s.append(gets())
    vc.append(getvc())
    fc.append(getfc())
    earned.append(earnings(p[-1],s[-1],vc[-1],fc[-1]))

ps = [s[i]*p[i] for i in range(n)]

fig,axs = plt.subplots(3)
fig.subplots_adjust(bottom=0.07,left=0.07)
axs[0].scatter( vc,earned,  linewidth=2.0)
axs[0].set_xlabel("Variable Costs")
axs[0].set_ylabel("Earnings")
axs[1].scatter( fc,earned,  linewidth=2.0)
axs[1].set_xlabel("Fixed Costs")
axs[1].set_ylabel("Earnings")
axs[2].scatter(ps ,earned,  linewidth=2.0)
axs[2].set_xlabel("Unit Sales x Unit Costs")
axs[2].set_ylabel("Earnings")


def lstsq(x,y):
    xbar = np.mean(x)
    ybar = np.mean(y)

    m = sum([(x[i]-xbar)*(y[i]-ybar) for i in range(len(x))])
    m/= sum((xi-xbar)**2 for xi in x)
    c = ybar-m*xbar
    return m,c

m1,c1 = lstsq(vc,earned)
m2,c2 = lstsq(fc,earned)
m3,c3 = lstsq(ps,earned)

axs[0].plot(vc,[c1+m1*x for x in vc],color= "red")
axs[1].plot(fc,[c2+m2*x for x in fc],color= "red")
axs[2].plot(ps,[c3+m3*x for x in ps],color= "red")

fig.show()

def sig(x):return sqrt(sum([(xi - np.mean(x))**2 for xi in x])/len(x))

Ssig1 = (sig(vc)*m1/sig(earned))**2
Ssig2 = (sig(fc)*m2/sig(earned))**2
Ssig3 = (sig(ps)*m3/sig(earned))**2

print(f"Sigma normalized derivative:\nvariable cost: {Ssig1}\nfixed cost: {Ssig2}\nprice x volume: {Ssig3}")
print(f"The sum of the normalized derivatives squated is: {Ssig1+Ssig2+Ssig3}")

#endregion

# wait for end of line
input("End of line. [ENTER]")