import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as spl
from scipy import integrate as int



#region Part1

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
print(bestAnt)
print(bestAnt.calcPathCost(mySim.dMatrix))


# for a in mySim.ants:
#     # print(a.path)
#     print(a.calcPathCost(mySim.dMatrix))

#endregion


#region Part2


#endregion
