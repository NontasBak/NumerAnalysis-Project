import matplotlib.pyplot as plt
import numpy as np


def LUmine(A):
    n = A.shape[0]
    L = np.matrix( np.identity(n) )
    U = A
    for j in range(0,n-1):
        for i in range(j+1,n):
            mult = A[i,j] / A[j,j]
            A[i, j+1:n] = A[i, j+1:n] - mult * A[j, j+1:n]
            U[i, j+1:n] = A[i, j+1:n]
            L[i,j] = mult
            U[i,j] = 0
    return L, U

def QRmine(A):
    n = A.shape[0]
    Q = np.matrix( np.zeros( (n,n) ) )
    for j in range(0, n):
        q = A[:,j]
        for i in range(0, j):
            length_of_leg = np.sum( A[:,j].T * Q[:,i])
            q = q - length_of_leg * Q[:,i]
        Q[:,j] = q / np.linalg.norm(q)
    R = Q.T * A
    return Q, R

'''
A = np.matrix([[1,1,0] , [1,0,1] , [0,1,1]])
r = np.linalg.qr(A)
print(r[0])
print(r[1])

test = QRmine(A)
print(test[0])
print(test[1])
'''

def Hilbert(n):
    H = np.matrix( np.zeros( (n,n) ) )
    for i in range(0,n):
        for j in range(0,n):
            H[i,j] = 1 / (i + j + 1)
    return H

'''
Hil = Hilbert(5)
print(Hil)
'''

def solveWithQR_Axb(A, b):
    QR = QRmine(A)
    x = np.linalg.inv(QR[1]) * QR[0].T * b
    return x

def solvewithLU_Axb(A, b):
    LU = LUmine(A)
    x = np.linalg.inv(LU[1]) * np.linalg.inv(LU[0]) * b
    return x

'''
n = 4
Hil = Hilbert(n)
b = np.ones((n, 1))
solution1 = solveWithQR_Axb(Hil, b)

print(solution1)
print(Hil * solution1)


solution2 = solvewithLU_Axb(Hil, b)
Hil = Hilbert(n) #H A=LU xalaei ton pinaka A

print(solution2)
print(Hil * solution2)
'''



