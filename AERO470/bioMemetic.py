import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as spl
from scipy import integrate as int



#region Part1

class Ant:
    alpha = 1
    beta = 2
    def __init__(self,node):
        self.node = node

    def __repr__(self):
        return f"Ant at node: {self.node}"
    
    def calcVMatrix(self,dMatrix):
        return np.array([[1/d if i != self.node and d!=0 else 0 for i,d in enumerate(row)] for row in dMatrix])
        # where d is an element of dMatrix, row is a row of dMatrix i is a colunm number
    
    def probOfUsingEdge(self,vMatrix,pMatrix):
        n = len(dMatrix)
        for i in range(n):
            t = vMatrix[self.node,i]# tau for this node t[r,s]
            # eta = 

pMatrix = np.ones(5) # pheramone matrix
dMatrix = np.array([[0,10,12,11,14],[10,0,13,15,8],[12,13,0,9,14],[11,15,9,0,16],[14,8,14,16,0]]) # distance matrix
myAnt = Ant(4)
print(myAnt.calcVMatrix(dMatrix))

    

#endregion


#region Part2


#endregion
