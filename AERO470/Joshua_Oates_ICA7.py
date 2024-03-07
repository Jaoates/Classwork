import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as spl
from scipy import integrate as int

A = np.array([[1,3,5],[2,5,1],[2,3,8]])
b = np.array([10,8,3])
print(f"A: \n {A}")
print(f"b: \n {b}")

Ar = np.linalg.matrix_rank(A)
print(f"rank(A): \n {Ar}")

try:
    Ai = np.invert(A)
except:
    Ai = "A is not invertible"
    
print(f"inv(A): \n {Ai}")

Ad = np.linalg.det(A)
print(f"det(A): \n {Ad}")

At = np.trace(A)
print(f"tr(A): \n {At}")

Ae = np.linalg.eig(A)
print(f"eig(A): \n {Ae}")

An = np.linalg.norm(A)
print(f"norm(A): \n {An}")

lu = spl.lu(A)
print(f"L,U: \n {lu}")

svd = np.linalg.svd(A)
print(f"L,U: \n {svd}")

x = np.linalg.solve(A,b)
print(f"x: {x}")

def vdp(t,X,mu):
    x = X[0]
    xd = X[1]
    xdd = mu*(1-x**2)*xd - x
    return np.array([xd,xdd])

mu = 1
tspan = [0,10]
t_eval=np.arange(tspan[0],tspan[1],.1)

X0 = []
xspan = np.arange(-3,3)
yspan = np.arange(-3,3)
for i in xspan:
    for j in yspan:
        X0.append( np.array([i,j]) )
ret = []


for i,X0i in enumerate(X0):
    ret.append(int.solve_ivp(lambda t,X: vdp(t,X,mu),tspan,X0i,t_eval=t_eval))
    

plt.style.use('_mpl-gallery')
fig,ax = plt.subplots()
margin = .1
fig.subplots_adjust(bottom=margin,left=margin,top = 1-margin,right = 1-margin)

for reti in ret:
    ax.plot(reti.y[0],reti.y[1])

# plt.title("Search space for BFA")
fig.show()


input("End of line [ENTER]")