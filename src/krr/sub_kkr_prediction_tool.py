import numpy as np
import os
import time

import sub_coord_single
import sub_kernel_regression_tool
from sklearn.metrics.pairwise import rbf_kernel
import sys
import shutil
from numba import jit

sys.dont_write_bytecode = True


class kkr_prediction_single():
    def __init__(self, n_x_dim, n_y_dim, x_train, para_kernel, para_gamma, para_alpha, fit_path, work_path, rescale,
                 label_grad):
        self.n_x_dim = n_x_dim
        self.n_y_dim = n_y_dim
        self.x_train = x_train
        self.para_kernel = para_kernel
        self.para_gamma = para_gamma
        self.para_alpha = para_alpha
        self.fit_path = fit_path
        self.work_path = work_path
        self.rescale = rescale
        self.label_grad = label_grad

    def kkr_prediction_energy_single(self):
        kkr = sub_kernel_regression_tool.kkr_all(self.n_x_dim, self.para_kernel, self.para_gamma, self.para_alpha,
                                                 self.fit_path, self.rescale)
#        print("Load all data.")
        kkr.load_data_train()
        kkr.load_data_test()
#        print("Scale data")
        kkr.rescale_data()
#        print("Test KKR")
        kkr.test_kkr()
#        print("Test error is: ", kkr.test_error)

    def kkr_prediction_energy_gradient(self):

        start_all = time.time()

        start = time.time()

        kkr = sub_kernel_regression_tool.kkr_all(self.n_x_dim, self.x_train, self.para_kernel, self.para_gamma,
                                                 self.para_alpha,
                                                 self.fit_path, self.rescale)
        elapsed = (time.time() - start)
        print("Initial kkr class Time used:",elapsed)

        start = time.time()

        kkr.load_data_train()

        elapsed = (time.time() - start)
        print("load_data_train Time used:",elapsed)

        start = time.time()

        kkr.load_data_prediction()

        elapsed = (time.time() - start)
        print("load_data_prediction Time used:",elapsed)

        start = time.time()

#        parameter_file = str(kkr.fit_path) + '/fitting_para.dat'
#        coeff = np.loadtxt(parameter_file)

        parameter_file = str(kkr.fit_path) + '/fitting_para.dat.npy'

        coeff = np.load(parameter_file)

#        print(coeff.shape)
        coeff = coeff.reshape(kkr.n_train, 1)
#        print(coeff.shape)

        kernel_term = rbf_kernel(kkr.x_pre, kkr.x_train, kkr.gamma)

        x_diff = np.arange(kkr.n_pre * kkr.n_train * kkr.n_x_dim).reshape(kkr.n_pre, kkr.n_train, kkr.n_x_dim)
        x_diff = x_diff.astype(float)

        y_pre_kkr = np.array([[0.0]])

        gradient_kkr = kkr_prediction_gradient(x_diff,kernel_term,coeff,kkr.n_train,kkr.x_train,kkr.x_pre,kkr.gamma,kkr.n_x_dim,kkr.n_pre,kkr.mean_y,kkr.scale_x_factor,kkr.x_pre_old,y_pre_kkr)

        filename = str(kkr.fit_path) + '/energy_gradient_col.dat'
        np.savetxt(filename, gradient_kkr)

        elapsed = (time.time() - start)
        print("prediction Time used:",elapsed)

        elapsed_all = (time.time() - start_all)
        print("Total predication Time Used",elapsed_all)


class single_x_descriptor():
    def __init__(self, filename, file_type, label_x_descriptor):
        self.file_name = filename
        self.file_type = file_type
        self.label_x_descriptor = label_x_descriptor

    def descriptor(self):
        n_atom, atom_name, coord = sub_coord_single.distance_matrix(self.file_name, self.file_type)
        self.n_atom = n_atom

    def extract_x_input(self, file_distance):

        x_all = np.loadtxt(file_distance)
        n_dim = x_all.shape[0]
        n_input_len = int((n_dim * (n_dim - 1)) / 2)
        x_input = np.zeros((1, n_input_len))

        i_input = 0
        for i_row in range(n_dim):
            for i_column in range(i_row):
                x_input[0, i_input] = x_all[i_row, i_column]
                i_input = i_input + 1
        self.x_input = x_input

    def descriptor_x(self):
        file_distance_matrix = 'distance_matrix.dat'
        self.extract_x_input(file_distance_matrix)
        x_old = self.x_input
        np.savetxt('x_input.dat', x_old)

    def descriptor_1_x(self):
        file_distance_matrix = 'distance_matrix.dat'
        self.extract_x_input(file_distance_matrix)
        x_old = 1.0 / self.x_input
        np.savetxt('x_inversed_input.dat', x_old)

    def descriptor_coulomb_matrix(self):
        file_distance_matrix = 'coulomb_matrix.dat'
        self.extract_x_input(file_distance_matrix)
        x_old = self.x_input
        np.savetxt('coulomb_input.dat', x_old)

    def descriptor_file(self):
        if self.label_x_descriptor == 0:
            filename = 'coulomb_input.dat'
        if self.label_x_descriptor == 1:
            filename = 'x_input.dat'
        if self.label_x_descriptor == 2:
            filename = 'x_inversed_input.dat'
        self.filename_pre = filename

    def xyz_to_descriptor(self):

        file_list = ['x_input.dat', 'x_inversed_input.dat', 'coulomb_matrix.dat', 'coulomb_input.dat', 'x_input.dat',
                     'x_inversed_input.dat', './x_pre.dat']
        for file_tmp in file_list:
            if os.path.exists(file_tmp):
                os.remove(file_tmp)

        self.descriptor()
        self.descriptor_x()
        self.descriptor_1_x()
        self.descriptor_coulomb_matrix()
        self.descriptor_file()

        curr_dir = os.getcwd()

        original_path = curr_dir + '/' + str(self.filename_pre)
        dest_path = curr_dir + '/' + 'x_pre.dat'

        shutil.copy2(original_path, dest_path)


class gradient_descriptor_to_xyz():
    def __init__(self, n_atom, n_x_dim, n_y_dim, label_x_descriptor):
        self.n_atom = n_atom
        self.n_x_dim = n_x_dim
        self.n_y_dim = n_y_dim
        self.label_x_descriptor = label_x_descriptor
        if self.label_x_descriptor != 0:
            print("The gradient calculations only support the coulumb matrix as input vector! ")
            raise IOError

    def read_gra_file(self):
        nx_dim = self.n_x_dim
        filename_col = './energy_gradient_col.dat'
        xx = np.loadtxt(filename_col)
        nn_xx_dim = xx.ndim
        if nn_xx_dim == 1:
            xx = xx.reshape(1, -1)
#        print(xx.ndim, xx.shape)
        self.qq = xx[0, 0:nx_dim]
        self.energy = xx[0, nx_dim]
        self.energy_pre = xx[0, nx_dim + 1]
        self.gradient_coul = xx[0, nx_dim + 2: 2 * nx_dim + 2]
        self.qq = self.qq.reshape(1, nx_dim)
        #        self.energy = self.energy.reshape(1, 1)
        self.gradient_coul = self.gradient_coul.reshape(1, nx_dim)

        n_dim = (self.n_atom * (self.n_atom - 1)) / 2
        if n_dim != nx_dim:
            print("The atomic number and the dimension of C matrix is not consistent!")
            raise IOError

    def read_geom_file(self):
        file_geom = "./standard.xyz"
        filein1 = open(file_geom, 'r')
        file_all = filein1.read()
        filein1.close()
        atom_name = []
        coord = np.zeros((self.n_atom, 3))
        charge = np.zeros((self.n_atom))
        file_text = file_all.split('\n')
        for i in range(self.n_atom):
            atom_name.append(file_text[i + 2].split()[0])
            coord[i, 0] = file_text[i + 2].split()[1]
            coord[i, 1] = file_text[i + 2].split()[2]
            coord[i, 2] = file_text[i + 2].split()[3]
            charge[i] = file_text[i + 2].split()[4]
        self.atom_name = atom_name
        self.coord = coord
        self.charge = charge

    def xyz_to_distance(self):
        n_atom = self.n_atom
        dis_2 = np.zeros((n_atom, n_atom))
        for i in range(n_atom):
            for j in range(n_atom):
                for k in range(3):
                    dis_2[i, j] = dis_2[i, j] + (self.coord[i, k] - self.coord[j, k]) ** 2
        self.dis_2 = dis_2

    def gra_coulumb_to_xyz(self):

        n_atom = self.n_atom
        gra_coul_2d = np.zeros((self.n_atom, self.n_atom))
        gra_coul_2d = gra_coul_2d.astype(float)
        qq_coul_2d = np.zeros((self.n_atom, self.n_atom))
        qq_coul_2d = qq_coul_2d.astype(float)

        i_input = 0
        for i_row in range(self.n_atom):
            for i_column in range(i_row):
                qq_coul_2d[i_row, i_column] = self.qq[0, i_input]
                qq_coul_2d[i_column, i_row] = qq_coul_2d[i_row, i_column]
                gra_coul_2d[i_row, i_column] = self.gradient_coul[0, i_input]
                gra_coul_2d[i_column, i_row] = gra_coul_2d[i_row, i_column]
                i_input = i_input + 1

        gradient_car = np.zeros((n_atom, 3))
        for i_atom in range(n_atom):
            for j_dim in range(3):
                for k_atom in range(n_atom):
                    if k_atom != i_atom:
                        gradient_car[i_atom, j_dim] = gradient_car[i_atom, j_dim] \
                                                      + gra_coul_2d[i_atom, k_atom] * \
                                                        qq_coul_2d[i_atom, k_atom] * \
                                                        1.0 / self.dis_2[i_atom, k_atom] * \
                                                        (self.coord[i_atom, j_dim] - self.coord[k_atom, j_dim])
        self.gradient_car = -  gradient_car

#    def energy_gradient_unit(self):
#        print('be careful for the unit conversion')
        #self.energy = self.energy / 27.2112
        #self.gradient_car = self.gradient_car / 27.2112
        #self.energy_pre = self.energy_pre / 27.2112

    def save_gradient(self):

        file_grad = './gradient.xyz'
        fileout1 = open(file_grad, 'w')
        fileout1.write(' ' + str(self.n_atom) + '\n')
        fileout1.write('Coordinate atom unit' + '\n')
        for i in range(self.n_atom):
            fileout1.write(str(self.atom_name[i]) + '   ' + "%14.8f" % (self.coord[i, 0]) + '   ' + "%14.8f" % (
                self.coord[i, 1]) + '   ' + "%14.8f" % (self.coord[i, 2]) + '   ' + str(self.charge[i]) + '\n')
        fileout1.write(' Energy of State:  ' + str(self.energy) + ' Hartree   ' + str(self.energy * 27.2112) + ' eV \n')
        fileout1.write(
            ' Energy of State:  ' + str(self.energy_pre) + ' Hartree   ' + str(self.energy_pre * 27.2112) + ' eV \n')
        fileout1.write(' Gradient atomic unit ' + '\n')
        for i in range(self.n_atom):
            fileout1.write('   ' + "%14.8f" % (self.gradient_car[i, 0]) + '   ' + "%14.8f" % (
                self.gradient_car[i, 1]) + '   ' + "%14.8f" % (self.gradient_car[i, 2]) + '  \n')
        fileout1.close()

        file_grad_1d = './gradient_1d.xyz'
        np.savetxt(file_grad_1d, self.gradient_car, fmt='%14.8f   ' * self.gradient_car.shape[1])

    def gradient_xyz(self):
        self.read_gra_file()
        self.read_geom_file()
        self.xyz_to_distance()
        self.gra_coulumb_to_xyz()
#        self.energy_gradient_unit()
        self.save_gradient()


class kkr_single_all_step():
    def __init__(self, n_x_dim, n_y_dim, x_train, para_kernel, para_gamma, para_alpha, fit_path, work_path, rescale,
                 label_grad,
                 label_x_descriptor):
        self.n_x_dim = n_x_dim
        self.n_y_dim = n_y_dim
        self.x_train = x_train
        self.para_kernel = para_kernel
        self.para_gamma = para_gamma
        self.para_alpha = para_alpha
        self.fit_path = fit_path
        self.work_path = work_path
        self.rescale = rescale
        self.label_grad = label_grad
        self.label_x_descriptor = label_x_descriptor

    def kkr_all_step(self):
        filename = 'geom.xyz'
        filetype = 'xyz'
        xyz_to_input_x = single_x_descriptor(filename, filetype, self.label_x_descriptor)
        xyz_to_input_x.xyz_to_descriptor()
        n_atom = xyz_to_input_x.n_atom

        kkr_pre = kkr_prediction_single(self.n_x_dim, self.n_y_dim, self.x_train, self.para_kernel, self.para_gamma,
                                        self.para_alpha,
                                        self.fit_path, self.work_path, self.rescale, self.label_grad)

        start = time.time()

        kkr_pre.kkr_prediction_energy_gradient()

        elapsed = (time.time() - start)
        print("Prediction of energies and gradients time used:",elapsed)

        convert = gradient_descriptor_to_xyz(n_atom, self.n_x_dim, self.n_y_dim, self.label_x_descriptor)
        convert.gradient_xyz()

    def kkr_all_step_numerical_gra(self):

        filetype = 'xyz'
        n_atom, atom_name, coord_ref = sub_coord_single.distance_matrix(file_name_zero, file_type)

        step_size = 0.05
        for i_atom in range(n_atom):
            for j in range(3):
                dir_tmp = str(i_atom + 1) + '_' + str(j + 1)
                os.mkdir(dir_tmp)
                coord = coord_ref
                coord[i_atom][j] = coord_ref[i_atom][j] + step_size

        xyz_to_input_x = single_x_descriptor(filename, filetype, self.label_x_descriptor)
        xyz_to_input_x.xyz_to_descriptor()

        kkr_pre = kkr_prediction_single(self.n_x_dim, self.n_y_dim, self.para_kernel, self.para_gamma, self.para_alpha,
                                        self.fit_path, self.work_path, self.rescale, self.label_grad)
        kkr_pre.kkr_prediction_energy_gradient()

@jit(nopython=True)
def kkr_prediction_gradient(x_diff,kernel_term,coeff,n_train,x_train,x_pre,gamma,n_x_dim,n_pre,mean_y,scale_x_factor,x_pre_old,y_pre_kkr):

#    kernel_term = rbf_kernel(x_pre, x_train, gamma)
#    print(kernel_term.shape)

#    x_diff = np.arange(n_pre * n_train * n_x_dim).reshape(n_pre, n_train, n_x_dim)
#    x_diff = x_diff.astype(float)
#    x_diff = float(x_diff)
#    print(x_diff.shape)

#    print("Compute energy")
    energy = np.dot(kernel_term, coeff)
    for i_pre in range(n_pre):
        energy[i_pre, :] = energy[i_pre, :] + mean_y[:]

#    print("Compute distance matrix")

#    start = time.clock()

    for i_dim in range(n_x_dim):
        for i_train in range(n_train):
            for i_test in range(n_pre):
                x_diff[i_test, i_train, i_dim] = x_pre[i_test, i_dim] - x_train[i_train, i_dim]
#    print(x_diff.shape)


    for i_dim in range(n_x_dim):
        x_diff[:, :, i_dim] = -2 * gamma * x_diff[:, :, i_dim] / scale_x_factor[i_dim]

#    elapsed = (time.clock() - start)
#    print("compute distance matrix Time used:",elapsed)

#    print("Compute KKR gradient")

#    start = time.clock()

    kkr_gradient = np.zeros((n_pre, n_x_dim))
    for i_dim in range(n_x_dim):
        for i_test in range(n_pre):
            for i_train in range(n_train):
                kkr_gradient[i_test, i_dim] = kkr_gradient[i_test, i_dim] + coeff[i_train, 0] * x_diff[
                    i_test, i_train, i_dim] * kernel_term[i_test, i_train]


#    elapsed = (time.clock() - start)
#    print("compute KKR gradient Time used:",elapsed)

    # Check the fitting for training data
    gradient_kkr = np.hstack((x_pre_old, energy, y_pre_kkr))
    gradient_kkr = np.hstack((gradient_kkr, kkr_gradient))

    return (gradient_kkr)

if __name__ == "__main__":
    print('Nothing Done!')
