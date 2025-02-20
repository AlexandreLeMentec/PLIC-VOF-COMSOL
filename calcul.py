from math import *

phi_tr = 0
phi_cr = 0

x0 = 5
z0 = 0
nx = -sqrt(2)/2
nz = sqrt(2)/2
delta_x = 1
delta_z = 1
c = sqrt(delta_x**2 + delta_z**2)
phi = 0.5

nuo = x0 / delta_x
mx = (abs(nz)*delta_z) / (abs(nz)*delta_z + (abs(nx)*delta_x))
mz = (abs(nx)*delta_x) / (abs(nz)*delta_z + (abs(nx)*delta_x))

print("######################################################")

print("nu0 = ", nuo)
print("mx = ", mx)
print("mz = ", mz)

if mx <= mz : 
    if nx >=0 : 
        print("Cas mx <= mz & nx >=0  ")
        phi_tr = (mx/(2.*mz)*(nuo+1./3.)/(nuo+1./2.))
        phi_cr = (2.*mx/(3.*mz)*(nuo**3)/(nuo+1./2.)) 
    else :
        print("Cas mx <= mz & nx < 0  ")
        phi_tr = (mz/(2.*mx)*(nuo+mz/(3.*mx))/(nuo+1./2.))
        phi_cr = (2.*mx/(3.*mz)*(nuo**3)/(nuo+1./2.)) 
else : 
    if nx >=0 : 
        print("Cas mx > mz & nx >=0  ")
        phi_tr = (mx/(2.*mz)*(nuo+mz/(3*mx))/(nuo+1./2.))
        phi_cr = (2.*mx/(3.*mz)*(nuo**3)/(nuo+1./2.)) 
    else  :
        print("Cas mx > mz & nx < 0  ")
        phi_tr = (mz/(2.*mx)*(nuo+1-mz/(3*mx))/(nuo+1./2.))
        phi_cr = (2.*mx/(3.*mz)*((nuo+1)**3)/(nuo+1./2.))

print("phi_tr = ", phi_tr)
print("phi_cr = ", phi_cr)

