from ctypes import *
import numpy as np
import subprocess as sub

filepath = '/mnt/c/Users/33614/Documents/GitHub/PLIC-VOF-COMSOL/plic.so'
try:
    sub.run(["gcc", "-fPIC", "-shared", "-o", "plic.so", "plic_vof.c"])
    print("Compilation successful!")
except sub.CalledProcessError as e:
    print(f"Compilation failed: {e}")


C = 0.5
nx = np.sqrt(2)/2
nz = nx
dx = 1.0
dz = 1.0
x = 4.0
z = 0.0

func = CDLL(filepath)
print('############ Valeur calculée ################')
func.plic_cyl.argtypes = [c_float, c_float, c_float, c_float, c_float, c_float, c_float]
print("Si = ", func.plic_cyl(C,nx,nz,dx,dz,x,z))
print('')
print('##############################################')