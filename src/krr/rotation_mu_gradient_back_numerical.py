#!/usr/bin/env python3
import numpy as np
import os
import os.path
#from numba import jit

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

#print(x1,y1,z1,x2,y2,z2,x3,y3,z3)

def get_rotation_matrix(x1,y1,z1,x2,y2,z2,x3,y3,z3):

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

   d1 = (b2*c3-c2*b3)/tmp  
   d2 = (c2*a3-a2*c3)/tmp
   d3 = (a2*b3-b2*a3)/tmp
   
   e1 = (c1*b3-b1*c3)/tmp
   e2 = (a1*c3-c1*a3)/tmp
   e3 = (b1*a3-a1*b3)/tmp

   f1 = (b1*c2-c1*b2)/tmp
   f2 = (c1*a2-a1*c2)/tmp
   f3 = (a1*b2-b1*a2)/tmp

#   print (d1,d2,d3,e1,e2,e3,f1,f2,f3)

#   rotation_matrix_inverse = rotation_matrix.I

   rotation_inverse_array = np.array([[d1,e1,f1],[d2,e2,f2],[d3,e3,f3]])

   rotation_matrix_inverse = np.matrix(rotation_inverse_array)

   return (rotation_matrix, rotation_matrix_inverse)


rotation_matrix, rotation_matrix_inverse = get_rotation_matrix(x1,y1,z1,x2,y2,z2,x3,y3,z3)

length = 0.001

rotation_matrix_x1_plus,  rotation_matrix_inverse_x1_plus  = get_rotation_matrix(x1+length,y1,z1,x2,y2,z2,x3,y3,z3)
rotation_matrix_x1_minux, rotation_matrix_inverse_x1_minus = get_rotation_matrix(x1-length,y1,z1,x2,y2,z2,x3,y3,z3)

rotation_matrix_y1_plus,  rotation_matrix_inverse_y1_plus  = get_rotation_matrix(x1,y1+length,z1,x2,y2,z2,x3,y3,z3)
rotation_matrix_y1_minux, rotation_matrix_inverse_y1_minus = get_rotation_matrix(x1,y1-length,z1,x2,y2,z2,x3,y3,z3)

rotation_matrix_z1_plus,  rotation_matrix_inverse_z1_plus  = get_rotation_matrix(x1,y1,z1+length,x2,y2,z2,x3,y3,z3)
rotation_matrix_z1_minux, rotation_matrix_inverse_z1_minus = get_rotation_matrix(x1,y1,z1-length,x2,y2,z2,x3,y3,z3)

rotation_matrix_x2_plus,  rotation_matrix_inverse_x2_plus  = get_rotation_matrix(x1,y1,z1,x2+length,y2,z2,x3,y3,z3)
rotation_matrix_x2_minux, rotation_matrix_inverse_x2_minus = get_rotation_matrix(x1,y1,z1,x2-length,y2,z2,x3,y3,z3)

rotation_matrix_y2_plus,  rotation_matrix_inverse_y2_plus  = get_rotation_matrix(x1,y1,z1,x2,y2+length,z2,x3,y3,z3)
rotation_matrix_y2_minux, rotation_matrix_inverse_y2_minus = get_rotation_matrix(x1,y1,z1,x2,y2-length,z2,x3,y3,z3)

rotation_matrix_z2_plus,  rotation_matrix_inverse_z2_plus  = get_rotation_matrix(x1,y1,z1,x2,y2,z2+length,x3,y3,z3)
rotation_matrix_z2_minux, rotation_matrix_inverse_z2_minus = get_rotation_matrix(x1,y1,z1,x2,y2,z2-length,x3,y3,z3)

rotation_matrix_x3_plus,  rotation_matrix_inverse_x3_plus  = get_rotation_matrix(x1,y1,z1,x2,y2,z2,x3+length,y3,z3)
rotation_matrix_x3_minux, rotation_matrix_inverse_x3_minus = get_rotation_matrix(x1,y1,z1,x2,y2,z2,x3-length,y3,z3)

rotation_matrix_y3_plus,  rotation_matrix_inverse_y3_plus  = get_rotation_matrix(x1,y1,z1,x2,y2,z2,x3,y3+length,z3)
rotation_matrix_y3_minux, rotation_matrix_inverse_y3_minus = get_rotation_matrix(x1,y1,z1,x2,y2,z2,x3,y3-length,z3)

rotation_matrix_z3_plus,  rotation_matrix_inverse_z3_plus  = get_rotation_matrix(x1,y1,z1,x2,y2,z2,x3,y3,z3+length)
rotation_matrix_z3_minux, rotation_matrix_inverse_z3_minus = get_rotation_matrix(x1,y1,z1,x2,y2,z2,x3,y3,z3-length)


dx1 = (rotation_matrix_inverse_x1_plus - rotation_matrix_inverse_x1_minus) / (2*length)
dx2 = (rotation_matrix_inverse_x2_plus - rotation_matrix_inverse_x2_minus) / (2*length)
dx3 = (rotation_matrix_inverse_x3_plus - rotation_matrix_inverse_x3_minus) / (2*length)

dy1 = (rotation_matrix_inverse_y1_plus - rotation_matrix_inverse_y1_minus) / (2*length)
dy2 = (rotation_matrix_inverse_y2_plus - rotation_matrix_inverse_y2_minus) / (2*length)
dy3 = (rotation_matrix_inverse_y3_plus - rotation_matrix_inverse_y3_minus) / (2*length)

dz1 = (rotation_matrix_inverse_z1_plus - rotation_matrix_inverse_z1_minus) / (2*length)
dz2 = (rotation_matrix_inverse_z2_plus - rotation_matrix_inverse_z2_minus) / (2*length)
dz3 = (rotation_matrix_inverse_z3_plus - rotation_matrix_inverse_z3_minus) / (2*length)

#print("dx1")
#print(dx1)
#print("dx2")
#print(dx2)
#print("dx3")
#print(dx3)
#print("dy1")
#print(dy1)
#print("dy2")
#print(dy2)
#print("dy3")
#print(dy3)
#print("dz1")
#print(dz1)
#print("dz2")
#print(dz2)
#print("dz3")
#print(dz3)

mu_gradient_all_X = np.loadtxt('gradient_00_x.dat')
mu_gradient_all_Y = np.loadtxt('gradient_00_y.dat')
mu_gradient_all_Z = np.loadtxt('gradient_00_z.dat')

mu_all_data = np.loadtxt('mu_all_00.dat')
mu_all = np.array([mu_all_data[0],  mu_all_data[1], mu_all_data[2]])

for i in range(10):

    mu_gradient_rotated_x = np.array([mu_gradient_all_X[i,0],  mu_gradient_all_Y[i,0], mu_gradient_all_Z[i,0]])
    mu_gradient_rotated_y = np.array([mu_gradient_all_X[i,1],  mu_gradient_all_Y[i,1], mu_gradient_all_Z[i,1]])
    mu_gradient_rotated_z = np.array([mu_gradient_all_X[i,2],  mu_gradient_all_Y[i,2], mu_gradient_all_Z[i,2]])
    
    if (i==0):

        gradient_rotate_matrix_dx = dx1
        gradient_rotate_matrix_dy = dy1
        gradient_rotate_matrix_dz = dz1

    if (i==1):
        gradient_rotate_matrix_dx = dx2
        gradient_rotate_matrix_dy = dy2
        gradient_rotate_matrix_dz = dz2

    if (i==3):
        gradient_rotate_matrix_dx = dx3
        gradient_rotate_matrix_dy = dy3
        gradient_rotate_matrix_dz = dz3

    elif (i==2 or i>3):

        gradient_rotate_matrix_dx = np.matrix(np.zeros((3,3)))
        gradient_rotate_matrix_dy = np.matrix(np.zeros((3,3)))
        gradient_rotate_matrix_dz = np.matrix(np.zeros((3,3)))


    mu_rotate_x_back = np.dot(mu_gradient_rotated_x,rotation_matrix_inverse) + np.dot(mu_all,gradient_rotate_matrix_dx)
    mu_rotate_y_back = np.dot(mu_gradient_rotated_y,rotation_matrix_inverse) + np.dot(mu_all,gradient_rotate_matrix_dy)
    mu_rotate_z_back = np.dot(mu_gradient_rotated_z,rotation_matrix_inverse) + np.dot(mu_all,gradient_rotate_matrix_dz)

    
    mu_gradient_all_X[i,0] = mu_rotate_x_back[0,0]
    mu_gradient_all_Y[i,0] = mu_rotate_x_back[0,1]
    mu_gradient_all_Z[i,0] = mu_rotate_x_back[0,2]
    
    mu_gradient_all_X[i,1] = mu_rotate_y_back[0,0]
    mu_gradient_all_Y[i,1] = mu_rotate_y_back[0,1]
    mu_gradient_all_Z[i,1] = mu_rotate_y_back[0,2]
    
    mu_gradient_all_X[i,2] = mu_rotate_z_back[0,0]
    mu_gradient_all_Y[i,2] = mu_rotate_z_back[0,1]
    mu_gradient_all_Z[i,2] = mu_rotate_z_back[0,2]



np.savetxt('gradient_00_x.dat_back',mu_gradient_all_X)
np.savetxt('gradient_00_y.dat_back',mu_gradient_all_Y)
np.savetxt('gradient_00_z.dat_back',mu_gradient_all_Z)

