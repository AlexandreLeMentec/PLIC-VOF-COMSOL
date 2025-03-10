import numpy as np

m = 5           # nz
n = 5           # nr
# Initialisation d'un tableau 10x10 rempli de zéros
A = np.zeros((m, n))
#init 
A[0][0] = 1
A[n-1][m-1] = n*m
for i in range(n-1):
    A[i+1][0] = A[i][0] + i+2
    A[n-i-2][m-1] = A[n-i-1][m-1] -(i+2)

for i in range(n-1) : 
    for j in range(1,m-i) : 
        A[i][j] = A[i][j-1] + (i+j)
        A[n-i-1][m-j-1] = A[n-i-1][m-j] -(i+j)

print(A)