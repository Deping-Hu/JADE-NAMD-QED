#!/usr/bin/env python
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

   rotation_inverse_array = np.array([[d1,e1,f1],[d2,e2,f2],[d3,e3,f3]])

   rotation_matrix_inverse = np.matrix(rotation_inverse_array)

   return (rotation_matrix, rotation_matrix_inverse)


rotation_matrix, rotation_matrix_inverse = get_rotation_matrix(x1,y1,z1,x2,y2,z2,x3,y3,z3)


###############
#### read the results before rotation
filename_in = 'ml_dipole.dat'
filename_out = 'ml_dipole_rotation_back.dat'

fpin = open(filename_in, "r")
fpout = open(filename_out, "w")

record = fpin.readline()
fpout.write(record)
record = fpin.readline()
fpout.write(record)

for i in range(10):
   record = fpin.readline()
   fpout.write(record)

record = fpin.readline()
fpout.write(record)

### get rotation mu
mu_x = []
mu_y = []
mu_z = []

for i in range(3):
   mu_x.append(float(fpin.readline()))
   mu_y.append(float(fpin.readline()))
   mu_z.append(float(fpin.readline()))

   mu_all = np.array([mu_x[i],  mu_y[i], mu_z[i]]) 

   mu_rotate_back = np.dot(mu_all,rotation_matrix_inverse) 

   fpout.write(str(mu_rotate_back[0,0])+'\n')
   fpout.write(str(mu_rotate_back[0,1])+'\n')
   fpout.write(str(mu_rotate_back[0,2])+'\n')

### get rotation mu gradient

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

length = length / 0.529

dx1 = (rotation_matrix_inverse_x1_plus - rotation_matrix_inverse_x1_minus) / (2*length)
dx2 = (rotation_matrix_inverse_x2_plus - rotation_matrix_inverse_x2_minus) / (2*length)
dx3 = (rotation_matrix_inverse_x3_plus - rotation_matrix_inverse_x3_minus) / (2*length)

dy1 = (rotation_matrix_inverse_y1_plus - rotation_matrix_inverse_y1_minus) / (2*length)
dy2 = (rotation_matrix_inverse_y2_plus - rotation_matrix_inverse_y2_minus) / (2*length)
dy3 = (rotation_matrix_inverse_y3_plus - rotation_matrix_inverse_y3_minus) / (2*length)

dz1 = (rotation_matrix_inverse_z1_plus - rotation_matrix_inverse_z1_minus) / (2*length)
dz2 = (rotation_matrix_inverse_z2_plus - rotation_matrix_inverse_z2_minus) / (2*length)
dz3 = (rotation_matrix_inverse_z3_plus - rotation_matrix_inverse_z3_minus) / (2*length)

fpin.readline()
fpout.write("Gradient of dipole\n")
for i in range(3):
    mu_all = np.array([mu_x[i],  mu_y[i], mu_z[i]]) 

    mu_gradient_all_X = []
    mu_gradient_all_Y = []
    mu_gradient_all_Z = []

    fpin.readline()
    for j in range(10):
       mu_gradient_all_X.append([])
       record = fpin.readline().split()
       mu_gradient_all_X[j].append(float(record[0]))
       mu_gradient_all_X[j].append(float(record[1]))
       mu_gradient_all_X[j].append(float(record[2]))

    fpin.readline()
    for j in range(10):
       mu_gradient_all_Y.append([])
       record = fpin.readline().split()
       mu_gradient_all_Y[j].append(float(record[0]))
       mu_gradient_all_Y[j].append(float(record[1]))
       mu_gradient_all_Y[j].append(float(record[2]))

    fpin.readline()
    for j in range(10):
       mu_gradient_all_Z.append([])
       record = fpin.readline().split()
       mu_gradient_all_Z[j].append(float(record[0]))
       mu_gradient_all_Z[j].append(float(record[1]))
       mu_gradient_all_Z[j].append(float(record[2]))


    for j in range(10):
    
        mu_gradient_rotated_x = np.array([mu_gradient_all_X[j][0],  mu_gradient_all_Y[j][0], mu_gradient_all_Z[j][0]])
        mu_gradient_rotated_y = np.array([mu_gradient_all_X[j][1],  mu_gradient_all_Y[j][1], mu_gradient_all_Z[j][1]])
        mu_gradient_rotated_z = np.array([mu_gradient_all_X[j][2],  mu_gradient_all_Y[j][2], mu_gradient_all_Z[j][2]])
        
        if (j==0):
    
            gradient_rotate_matrix_dx = dx1
            gradient_rotate_matrix_dy = dy1
            gradient_rotate_matrix_dz = dz1
    
        if (j==1):
            gradient_rotate_matrix_dx = dx2
            gradient_rotate_matrix_dy = dy2
            gradient_rotate_matrix_dz = dz2
    
        if (j==3):
            gradient_rotate_matrix_dx = dx3
            gradient_rotate_matrix_dy = dy3
            gradient_rotate_matrix_dz = dz3
    
        elif (j==2 or j>3):
    
            gradient_rotate_matrix_dx = np.matrix(np.zeros((3,3)))
            gradient_rotate_matrix_dy = np.matrix(np.zeros((3,3)))
            gradient_rotate_matrix_dz = np.matrix(np.zeros((3,3)))


        mu_rotate_x_back = np.dot(mu_gradient_rotated_x,rotation_matrix_inverse) + np.dot(mu_all,gradient_rotate_matrix_dx)
        mu_rotate_y_back = np.dot(mu_gradient_rotated_y,rotation_matrix_inverse) + np.dot(mu_all,gradient_rotate_matrix_dy)
        mu_rotate_z_back = np.dot(mu_gradient_rotated_z,rotation_matrix_inverse) + np.dot(mu_all,gradient_rotate_matrix_dz)

        mu_gradient_all_X[j][0] = mu_rotate_x_back[0,0]
        mu_gradient_all_Y[j][0] = mu_rotate_x_back[0,1]
        mu_gradient_all_Z[j][0] = mu_rotate_x_back[0,2]
        
        mu_gradient_all_X[j][1] = mu_rotate_y_back[0,0]
        mu_gradient_all_Y[j][1] = mu_rotate_y_back[0,1]
        mu_gradient_all_Z[j][1] = mu_rotate_y_back[0,2]
        
        mu_gradient_all_X[j][2] = mu_rotate_z_back[0,0]
        mu_gradient_all_Y[j][2] = mu_rotate_z_back[0,1]
        mu_gradient_all_Z[j][2] = mu_rotate_z_back[0,2]


    state_index = i*3+1

    fpout.write(" State: "+str(state_index)+"\n")
    for j in range(10):
        fpout.write('   ' + "%14.8f" % (mu_gradient_all_X[j][0]) + '   ' + "%14.8f" % (mu_gradient_all_X[j][1]) + '   ' + "%14.8f" % (mu_gradient_all_X[j][2]) + '  \n')
#
    state_index = i*3+2
    fpout.write(" State: "+str(state_index)+"\n")
    for j in range(10):
        fpout.write('   ' + "%14.8f" % (mu_gradient_all_Y[j][0]) + '   ' + "%14.8f" % (mu_gradient_all_Y[j][1]) + '   ' + "%14.8f" % (mu_gradient_all_Y[j][2]) + '  \n')
#
    state_index = i*3+3
    fpout.write(" State: "+str(state_index)+"\n")
    for j in range(10):
        fpout.write('   ' + "%14.8f" % (mu_gradient_all_Z[j][0]) + '   ' + "%14.8f" % (mu_gradient_all_Z[j][1]) + '   ' + "%14.8f" % (mu_gradient_all_Z[j][2]) + '  \n')
#
fpin.close()
fpout.close()

