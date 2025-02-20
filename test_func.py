from ctypes import *
import numpy as np
import subprocess as sub

filepath = '/mnt/c/Users/33614/Documents/GitHub/PLIC-VOF-COMSOL/plic.so'
try:
    sub.run(["gcc", "-fPIC", "-shared", "-o", "plic.so", "plic_vof.c"])
    print("Compilation successful!")
except sub.CalledProcessError as e:
    print(f"Compilation failed: {e}")

C = 0.95
nx = np.cos(7.6*np.pi/4.0)
nz = np.sin(7.6*np.pi/4.0)
dx = 1.0
dz = 1.0
x = 25.0
z = 0.0

func = CDLL(filepath)
print('############ Valeur calculée ################')
func.plic_cyl.argtypes = [c_float, c_float, c_float, c_float, c_float, c_float, c_float]
func.plic_cyl.restype = c_float
print("Si = ", func.plic_cyl(C,nx,nz,dx,dz,x,z))
print('')
print('##############################################')