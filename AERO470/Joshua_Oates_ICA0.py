# 1
x=5
print(x)

# 2
L = []
for i in range(5):
    L.append(i)
print(L)    

# 3
L2 = [a*x for a in L]
print(L2)

# 4
chars = ['c','b','a','d','e']
ords = [ord(c) for c in chars]

D = {c:o for (c,o) in zip(chars,ords)}
print(D)

D['!'] = 2
print(D)

# 5
D = {c:o for (c,o) in zip(sorted(D),[D[d] for d in sorted(D)])}
print(D)


