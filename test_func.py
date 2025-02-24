from ctypes import *
import numpy as np
import subprocess as sub
import matplotlib.pyplot as plt


C_mat_cas_1 = np.array([[1, 1, 0.6, 0.2, 0],
                        [1, 0.8, 0.1, 0, 0],
                        [1, 0.7, 0, 0, 0],
                        [1, 1, 0.5, 0, 0],
                        [1, 1, 1, 0.3, 0]])

C_mat_cas_2 = np.array([[0,0,0,0,0.5,1,1,1,1],
                        [0,0,0,0,0.5,1,1,1,1],
                        [0,0,0,0,0.5,1,1,1,1],
                        [0,0,0,0,0.25,0.5,0.5,0.5,0.5],
                        [0.5,0.5,0.25,0,0,0,0,0,0],
                        [1,1,0.5,0,0,0.2,0.2,0,0],
                        [1,1,0.5,0,0,0.2,0.2,0,0],
                        [1,1,0.75,0.25,0,0,0,0,0],
                        [1,1,1,0.5,0.3,0,0,0,0]]).T


def one_test(C, nx, nz, dx, dz, x, z):
    filepath = './plic.so'
    try:
        sub.run(["gcc", "-fPIC", "-shared", "-o", "plic.so", "plic_vof.c"])
        print("Compilation successful!")
    except sub.CalledProcessError as e:
        print(f"Compilation failed: {e}")
    func = CDLL(filepath)
    print('############ Valeur calculée ################')
    func.plic_cyl.argtypes = [c_float, c_float, c_float, c_float, c_float, c_float, c_float]
    func.plic_cyl.restype = c_float
    print("Si = ", func.plic_cyl(C,nx,nz,dx,dz,x,z))
    print('')
    print('##############################################')

def multi_test(C,dx,dz):
    filepath = './plic.so'
    try:
        sub.run(["gcc", "-fPIC", "-shared", "-o", "plic.so", "plic_vof.c"])
        print("> Compilation successful!")
    except sub.CalledProcessError as e:
        print(f"> Compilation failed: {e}")
    func = CDLL(filepath)
    print('> Calcul imputs ')
    func.plic_cyl.argtypes = [c_float, c_float, c_float, c_float, c_float, c_float, c_float]
    func.plic_cyl.restype = c_float
    x = np.zeros(np.shape(C))
    print(np.shape(C))
    z = np.zeros(np.shape(C))
    S = np.zeros(np.shape(C))
    grad = np.sqrt(np.pow(np.gradient(C,dx)[0],2) + np.pow(np.gradient(C,dz)[1],2))
    nx = -np.gradient(C,dx)[0]/np.abs(grad)
    nz = -np.gradient(C,dz)[1]/np.abs(grad)
    for i in range(np.shape(C)[0]):
        for j in range(np.shape(C)[1]):
            x[i, j] = dx/2 + dx*i
            z[i, j] = dz/2 + dz*j
    print('> Calcul imputs')
    print(' ')
    print('#################### CALCUL S: ###################')
    print(' ')
    for i in range(np.shape(C)[0]):
        for j in range(np.shape(C)[1]):
            print('Cell (',i,',',j,')----------------------------------')
            print('C = ', C[i, j])
            S[i,j] = func.plic_cyl(C[i,j],nx[i,j],nz[i,j],dx,dz,x[i,j],z[i,j])
            print('S =', S[i,j])
            if(C[i,j] != 0 and C[i,j]!= 1):
                if(S[i,j] == 0):
                    print(' PROBLEME #############################################################################')
            print('----------------------------------------------------')
    print('#################### Plot: ###################')
    # Create meshgrid for the x and y values (if they are 1D)
    

    # Plotting the heatmap
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(x, z, S, cmap='viridis')  # You can change the number of contour levels and colormap
    plt.colorbar(label='S values')  # Colorbar to show value scale
    plt.xlabel('r-axis')
    plt.ylabel('Z-axis')
    plt.title('Heatmap of Surface')
    plt.show()

multi_test(C_mat_cas_2, 1.0, 1.0)


