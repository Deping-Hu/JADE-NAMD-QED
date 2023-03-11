#!/usr/bin/env python
import numpy as np
import os
import os.path

filename = 'geom.xyz'

fpin = open(filename, "r")

# number of atom, 
line = fpin.readline()
natom = int(line)
line = fpin.readline()[0:-1]
molinfo = {'n_atom': natom, 'title': line}

atom = []
for j in range(natom):
    line = fpin.readline()
    rec = line.split()
    atomname, x, y, z = rec[0:4]
    record = {'name': atomname, 'coord': [float(x), float(y), float(z)]}
    atom.append(record)
mol = {'info': molinfo, 'atoms': atom}

fpin.close()

#print(mol['atoms'][0])

coord_all = mol['atoms']

x1 = coord_all[0]['coord'][0] 
y1 = coord_all[0]['coord'][1] 
z1 = coord_all[0]['coord'][2] 

x2 = coord_all[1]['coord'][0] 
y2 = coord_all[1]['coord'][1] 
z2 = coord_all[1]['coord'][2] 

x3 = coord_all[3]['coord'][0] 
y3 = coord_all[3]['coord'][1] 
z3 = coord_all[3]['coord'][2] 


tmp = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2

tmp = tmp**0.5

a1 =  x2-x1
a1 = a1/tmp

b1 =  y2-y1
b1 = b1/tmp

c1 =  z2-z1
c1 = c1/tmp

m = y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2)
n = z1*(x2-x3) + z2*(x3-x1) + z3*(x1-x2)
k = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
l = -x1*(y2*z3-y3*z2) -x2*(y3*z1-y1*z3) -x3*(y1*z2-y2*z1)

tmp = m**2 + n**2 + k**2
tmp = tmp**0.5

a3 = m/tmp
b3 = n/tmp
c3 = k/tmp

a2 = b1*c3-b3*c1
b2 = a3*c1-a1*c3
c2 = a1*b3-a3*b1

rotation_array = np.array([[a1,b1,c1],[a2,b2,c2],[a3,b3,c3]])
rotation_matrix = np.matrix(rotation_array)

tmp = a1*(b2*c3-c2*b3) - a2*(b1*c3-c1*b3) + a3*(b1*c2-c1*b2)

d1 = b2*c3-c2*b3
d2 = c2*a3-a2*c3
d3 = a2*b3-b2*a3

e1 = c1*b3-b1*c3
e2 = a1*c3-c1*a3
e3 = b1*a3-a1*b3

f1 = b1*c2-c1*b2
f2 = c1*a2-a1*c2
f3 = a1*b2-b1*a2

#   rotation_matrix_inverse = rotation_matrix.I

rotation_inverse_array = np.array([[d1,e1,f1],[d2,e2,f2],[d3,e3,f3]])
rotation_inverse_array = rotation_inverse_array/tmp
rotation_matrix_inverse = np.matrix(rotation_inverse_array)


mu_x_all = np.loadtxt('kkr_mu_x.dat')
mu_y_all = np.loadtxt('kkr_mu_y.dat')
mu_z_all = np.loadtxt('kkr_mu_z.dat')

for i in range(4):

    mu_rotated = np.array([mu_x_all[i],  mu_y_all[i], mu_z_all[i]])

#    print('mu_rotated')
#    print(mu_rotated)

    mu_rotate_back = np.dot(mu_rotated,rotation_matrix_inverse)

#    print('mu_rotate_back')
#    print(mu_rotate_back)

    mu_x_all[i] = mu_rotate_back[0,0]
    mu_y_all[i] = mu_rotate_back[0,1]
    mu_z_all[i] = mu_rotate_back[0,2]

np.savetxt('kkr_mu_x.dat_rotate_back',[mu_x_all])
np.savetxt('kkr_mu_y.dat_rotate_back',[mu_y_all])
np.savetxt('kkr_mu_z.dat_rotate_back',[mu_z_all])


