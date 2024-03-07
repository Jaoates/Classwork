import matplotlib.pyplot as plt
import numpy as np
import math as m




#region part1

values = np.array([23, 21, 8, 1, 3, 7, 18, 19, 17, 15, 24, 22, 6, 28, 4, 2, 27, 20, 5, 10])
weights = np.array([7, 2, 6, 9, 1, 5, 6, 1, 3, 4, 7, 9, 3, 7, 3, 4, 5, 1, 5, 4]) 
maxWeight = 45
n = 20
N = 2**20

print(f"for a selection of {n} items, there are {N} combinations in binary thats {bin(N-1)}")

def binlist(combo):
    myString = bin(combo)[2:].rjust(20,'0')
    # print(bin(combo).rjust(20,'0'))
    # print(myString)
    return [int(s) for s in myString]
# [binlist(b) for b in range(10)]

# combos = np.arange(N)
# def calcCombo(combo,numbers):
#     myString = bin(combo)[2:-0].rjust(20,'0')
#     return sum([int(myString[i])*numbers[i] for i in range(n)])

def calcCombo(combo,numbers):
    myList = binlist(combo)
    return sum([myList[i]*numbers[i] for i in range(n)])

def calcInt(blist):
    return sum([2**i * blist[i] for i in range(len(blist))])


# def calcCombo(combo,numbers):
#     return sum([(numbers[i]*((combo^(1<<i))>>i)) for i in range(len(numbers))])

v = []
W = []
I = []
# w = []

for i in range(N):
    # c = combos[i]
    w=(calcCombo(i,weights))
    # v.append(calcCombo(i,values))
    if w <= maxWeight:
        v.append(calcCombo(i,values))
        I.append(i)
        W.append(w)
    if (i%(N//10)) == 0:
        print(f"{100*i/N:3.1f}%")

print(len(v))
V = max(v)
i=v.index(V)
i = I[i]
w = calcCombo(i,weights)
# print(f"{v}\n{I}")
print(f"The maximum possible value is ${V} and is comprised of the combo {bin(i)}, the weight of this combo is {w} kg")


def showCombo(combo):
    b = binlist(combo)
    for i in range(len(b)):
        print(f"{b[i]}  |  {weights[i]}  |  {values[i]}")
    print("----------------------")
    print(f"   |  {w} |  {V}")


showCombo(i)

plt.style.use('_mpl-gallery')
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.06,left=0.07,top = 0.90)
ax.scatter(W,v,s = 1)
plt.title("Search space for BFA")
fig.show()

#endregion --------------------------------------------------------

#region Part 2
npts = 10000
searchSpace = [np.random.randint(0,N) for i in range(npts)]
v = []
W = []
I = []
w = float('nan')

for i in searchSpace:
    # c = combos[i]
    w=(calcCombo(i,weights))
    # v.append(calcCombo(i,values))
    if w <= maxWeight:
        v.append(calcCombo(i,values))
        I.append(i)
        W.append(w)

print(len(v))
V = max(v)
i=v.index(V)
i = I[i]
w = calcCombo(i,weights)
# print(f"{v}\n{I}")
print(f"The maximum possible value in the space searched is ${V} and is comprised of the combo {bin(i)}, the weight of this combo is {w} kg")

showCombo(i)

plt.style.use('_mpl-gallery')
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.06,left=0.07,top = 0.90)
ax.scatter(W,v,s = 1)
plt.title("Search space for MBFA")
fig.show()


plt.style.use('_mpl-gallery')
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.06,left=0.07,top = 0.90)
ax.plot(v)
plt.title("Solution fitness of each iteration MBFA")
fig.show()
#endregion --------------------------------------------------------

#region Part 3

genSize = 60
numGens = 10
pop = [np.random.randint(0,N) for i in range(genSize)]
mutationRate = 1
numCrossovers = 3

def mutate(combo,mutRate):
    for i in range(mutRate):
        combo = combo ^ (1<<np.random.randint(0,20))
    return combo

def crossover(combos):
    combo1 = np.random.choice(combos)
    combo2 = np.random.choice(combos)
    b1 = binlist(combo1)
    b2 = binlist(combo2)
    i = np.random.randint(0,len(b1))
    # print(i)
    # print(b1[0:i])
    # print(b2[i:])
    # print(b1[0:i]+b2[i:])
    return calcInt(b1[0:i]+b2[i:])

V = []
W = []
fatCombos =[]
portionToTake = .1
# bestCombo = 0
topSol = []
for i in range(numGens):
    print()
    combos = []
    for j in range(genSize):
        w = (calcCombo(pop[j],weights))
        if w<=maxWeight:
            combos.append(pop[j])
            fatCombos.append(pop[j])
        V = [calcCombo(c,values) for c in combos]
    # print(f"V first is {V}")
    v = V.copy()
    v.sort()
    v = v[-int(np.floor(portionToTake*genSize)):]
    print(f"V sorted is {v}")
    # print(f"best combo is {bestCombo}")
    # print(calcCombo(bestCombo,values))
    i = [V.index(vi) for vi in v]
    thinCombos = [combos[ii] for ii in i]
    topSol.append(v[-1])
    # if calcCombo(bestCombo,values) < v[-1]:
    #     bestCombo = thinCombos[-1]
    combos = np.random.choice(thinCombos,genSize)
    for j in range(numCrossovers):
        combos[np.random.randint(0,len(combos))] = crossover(combos)
    pop = [mutate(c,mutationRate) for c in combos]
    # print(pop)
    # pop.append(bestCombo)
    # print(pop)



v = [calcCombo(c,values) for c in fatCombos]
w = [calcCombo(c,weights) for c in fatCombos]
plt.style.use('_mpl-gallery')
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.06,left=0.07,top = 0.90)
ax.scatter(w,v,s = 1)
plt.title("Search space of GA")
fig.show()

plt.style.use('_mpl-gallery')
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.06,left=0.07,top = 0.90)
ax.plot(topSol)
plt.title("Solution fitness of each iteration GA")
fig.show()

#endregion --------------------------------------------------------
input("End of line [ENTER]")