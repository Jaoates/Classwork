import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as spl
from scipy import integrate as int
import copy as copy

regionsToRun = [1,2]


#region Part1
if 1 in regionsToRun:
    class Ant:
        alpha = 1
        beta = 2
        Q = 1
        nextNode = None
        def __init__(self,node: int) -> None:
            self.node = node
            self.path = [self.node]

        def __repr__(self):
            return f"Ant {self.path}"
        
        def calcVMatrix(self,dMatrix: np.ndarray) -> np.ndarray:
            return np.array([[1/d if i not in self.path and d!=0 else 0 for i,d in enumerate(row)] for row in dMatrix])
            # where d is an element of dMatrix, row is a row of dMatrix i is a colunm number
        
        def probOfNextNode(self,dMatrix: np.ndarray,pMatrix: np.ndarray) -> np.ndarray:
            n = len(dMatrix)
            assert(pMatrix.shape == (n,n))
            assert(dMatrix.shape == (n,n))
            p = np.zeros(n)
            vMatrix = self.calcVMatrix(dMatrix)
            for i in range(n):
                eta = vMatrix[self.node,i]# tau for this node t[r,s]
                tau = pMatrix[self.node,i]
                p[i] = tau**self.alpha * eta**self.beta
            p /= np.sum(p)
            return p

        def pickNextNode(self,dMatrix: np.ndarray,pMatrix: np.ndarray) -> None:
            n = len(dMatrix)
            assert(pMatrix.shape == (n,n))
            assert(dMatrix.shape == (n,n))
            if len(self.path) < n:
                p = self.probOfNextNode(dMatrix,pMatrix)
                p = np.cumsum(p)
                r = np.random.rand()
                for i,pi in enumerate(p):
                    if r <= pi:
                        self.nextNode = i
                        return None
                self.nextNode = n-1 # this line is cuaseing headaches, I dont think it should ever actually run, but it does?
            else:
                self.nextNode = self.path[0]


        def dropPheromone(self,dMatrix: np.ndarray,pMatrix: np.ndarray) -> np.ndarray:
            pMatrix[self.node,self.nextNode] += self.Q/dMatrix[self.node,self.nextNode]

        def updatePath(self):
            self.path.append(self.nextNode)
            self.node = self.nextNode
            self.nextNode = None

        def calcPathCost(self,dMatrix):
            n = len(self.path)
            return sum([dMatrix[self.path[i],self.path[i+1]] for i in range(n-1)])

    class AcoSim:    
        def __init__(self,nAnts: int, nRuns: int, rho ,dMatrix: np.ndarray,pMatrix: np.ndarray) ->None:
            self.pMatrix = pMatrix
            self.dMatrix = dMatrix
            self.nRuns = nRuns
            self.nAnts = nAnts
            self.rho = rho
            self.n = len(dMatrix)

        def spawnAnts(self):
            # self.ants = [Ant(np.random.randint(0,self.n)) for i in range(self.nAnts)]
            self.ants = [Ant(0) for i in range(self.nAnts)]

        def __repr__(self):
            return f"AcoSim with {len(self.pMatrix)} nodes"
        
        def evapPheromone(self):
            self.pMatrix = (1-self.rho)*self.pMatrix

        def rankAnts(self):
            self.ants = sorted(self.ants,key = lambda x: x.calcPathCost(self.dMatrix))
        
        def doSim(self):
            for i in range(self.nRuns):
                self.spawnAnts()
                for j in range(self.n): # ant should only step to the number of cities
                    for a in self.ants:
                        a.pickNextNode(self.dMatrix,self.pMatrix)
                        a.dropPheromone(self.dMatrix,self.pMatrix)
                        a.updatePath()
                        self.evapPheromone()
            self.rankAnts()
            bestAnt = self.ants[0]
            return bestAnt
                # print(self.ants)
                # print(self.pMatrix)        
                
    pMatrix = np.ones([5,5]) # pheromone matrix
    dMatrix = np.array([[0,10,12,11,14],
                        [10,0,13,15,8],
                        [12,13,0,9,14],
                        [11,15,9,0,16],
                        [14,8,14,16,0]]) # distance matrix

    mySim=AcoSim(3,20,.005,dMatrix,pMatrix)
    bestAnt = mySim.doSim()
    print("the best path is:")
    print(bestAnt)
    print("the value is:")
    print(bestAnt.calcPathCost(mySim.dMatrix))

#endregion


#region Part2
    
def ackleys(x,**kwargs):
    # ackleys function
    # from: https://www.sfu.ca/~ssurjano/ackley.html
    a = kwargs.get('a',20)
    b = kwargs.get('b',.2)
    c = kwargs.get('c',2*np.pi)
    d = len(x)
    return -a*np.exp(-b*np.sqrt((1/d)*sum([x**2 for x in x])))  -np.exp((1/d)*sum([np.cos(c*x) for x in x]))  +a +np.exp(1) 


class Particle:
    # note this class is highly inspired by Dr. Meheil's renowned code
    weight = .04 # inertia of particle
    cognitive = .5 # cognitive weight
    social = .6 # social weight
    def __init__(self,psosim) -> None:
        self.dom = psosim.domain
        self.dim = psosim.dimention
        self.pos = ( np.random.rand(self.dim) - .5 ) * 2 * self.dom
        self.bestPos = copy.deepcopy(self.pos) # im concerned this will create an alias... need deepCopy?
        self.vel = ( np.random.rand(self.dim) - .5 ) * 2 * self.dom 
        self.updateBest() # this will init a randomly generated particle with the correct value for its pos

    def __repr__(self):
        return f"Particle    val:{self.val}    pos:{self.pos}"

    def updatePos(self):
        self.pos += self.vel
        for i in range(self.dim):
            if self.pos[i] > self.dom:
                self.pos[i] = self.dom

    def updateVel(self,psosim):
        self.vel *= self.weight
        self.vel += self.cognitive * (self.bestPos - self.pos) 
        self.vel += self.social * (psosim.bestParticle.bestPos - self.pos)

    def updateBest(self):
        self.val = ackleys(self.pos)
        if self.val < ackleys(self.bestPos):
            self.bestPos = self.pos

class PsoSim():
    def __init__(self,nParticles,nRuns,simDomain,simDimention):
        self.nParticles = nParticles
        self.nRuns = nRuns
        self.domain = simDomain # what is the min and max distance in all dims of the sim
        self.dimention = simDimention # how many dimensions do the particles live in
        self.particles = [Particle(self) for i in range(self.nParticles)]
        self.updateBest()
    
    def updateBest(self):
        self.bestParticle = sorted(self.particles, key = lambda x: x.val)[0]

    def doSim(self):
        for i in range(self.nRuns):
            for p in self.particles:
                p.updateVel(self)
                p.updatePos()
                p.updateBest()
                self.updateBest()
        return(self.bestParticle)

sim = PsoSim(10,4,10,2)
bestParticle = sim.doSim()
print(f"The best particle is: {bestParticle}")

nx, ny = (300, 300)
x = np.linspace(-10, 10, nx)
y = np.linspace(-10, 10, ny)
X, Y = np.meshgrid(x, y)

Z = np.array([ackleys([x,y]) for (x,y) in zip(X.ravel(), Y.ravel())]).reshape(X.shape)
plt.contour(X,Y,Z)
plt.scatter([p.pos[0] for p in sim.particles],[p.pos[1] for p in sim.particles])
plt.show()

pass

#endregion
