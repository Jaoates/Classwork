
# 1
def myFactorial(n):
    f = 1
    n += 1
    for i in range(1,n):
        f*=i
    return f

print(myFactorial(1))
print(myFactorial(2))
print(myFactorial(3))
print(myFactorial(4))
print(myFactorial(0))
print(myFactorial(-1))
# print(myFactorial(4.3))

# 2

def matrixProduct(A,B):
    n = len(A) # num rows
    m = len(A[0]) # num cols
    l = [0]*n
    C = [[0]*n]*n
    c = [[0]*n]*n
    for i in range(n):
        for j in range(n):
            s = 0
            for k in range(m):
                s += A[j][k]*B[k][i]
            l[j] = s
        C[i] = l[:]
    for i in range(n):
        for j in range(n):
            l[j] = C[j][i]
        c[i] = l[:]
    return c

A = [[1,2,3],[4,5,6]]
B = [[10,11],[20,21],[30,31]]
print(A)
print(B)

C = matrixProduct(A,B)
print(C)

# 3
def elementWiseProduct(A,B):
    return [a*b for i, (a,b) in enumerate(zip(A,B))]

mass = [1.7, 4.2, 2.6, 5.4]
value = [2, 3, -1, 5]

print(elementWiseProduct(mass,value))