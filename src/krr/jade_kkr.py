#!/usr/bin/env python3 

import os

import sub_inp_json
import tools
import time
import numpy as np
import sys
import shutil
from sub_kkr_prediction_tool import kkr_single_all_step


sys.dont_write_bytecode = True

filename = 'kkr.input'
xxx = sub_inp_json.read_dat_with_label(filename)
sub_inp_json.dump_json('kkr.json', xxx)

ic = tools.interface_converter(filename='qm_interface')
mol = ic['mol']
natom = int(ic['mol']['natom'])
atoms = ic['mol']['atoms']
file_geom = open('geom.xyz', 'w')
file_geom.write('  ' + str(natom) + '\n')
file_geom.write('Ang\n')
for i_atom in range(natom):
    name = atoms[i_atom]['name']
    coord = atoms[i_atom]['coord']
    coord[0] = 0.529 * float(coord[0])
    coord[1] = 0.529 * float(coord[1])
    coord[2] = 0.529 * float(coord[2])
    file_geom.write(str(name) + '  ' + str(coord[0]) + '  ' + str(coord[1]) + '  ' + str(coord[2]) + '\n')
file_geom.close()


def kkr_single_geom_all():
    kkrname = "kkr.json"
    xxx_input = sub_inp_json.load_json(kkrname)
    n_x_dim = int(xxx_input['n_x_dim'])
    n_y_dim = int(xxx_input['n_y_dim'])
    label_x_descriptor = int(xxx_input['label_x_descriptor'])
    para_kernel = xxx_input['para_kernel']
    rescale = xxx_input['rescale']
    label_grad = xxx_input['label_grad']
    para_alpha = float(xxx_input['para_alpha'])
    para_gamma = float(xxx_input['para_gamma'])
    energy_zero = float(xxx_input['energy_zero'])

    curr_dir = os.getcwd()

    kkr_path = curr_dir + '/kkr'

    x_train = np.load('./kkr/x_train.dat.npy')

    for i_y_dim in range(n_y_dim):
        work_path = kkr_path + '/' + str(i_y_dim + 1)
        os.chdir(work_path)
        for i_file in ['coulomb_input.dat', 'distance_matrix.dat', 'gradient_1d.xyz', 'x_inversed_input.dat', 'coulomb_matrix.dat', 'energy_gradient_col.dat', 'geom.xyz', 'gradient.xyz', 'standard.xyz', 'x_input.dat', 'x_pre.dat']:

           file_path = work_path + '/' + i_file

           if os.path.isfile(file_path):

             os.remove(file_path)


        file_path = curr_dir + '/' + 'geom.xyz'

        shutil.copy2(file_path, work_path)

        fit_path = work_path

        kkr_single_all = kkr_single_all_step(n_x_dim, n_y_dim, x_train, para_kernel, para_gamma, para_alpha, fit_path,
                                             work_path,
                                             rescale, label_grad, label_x_descriptor)
        kkr_single_all.kkr_all_step()
        os.chdir(curr_dir)

    file_kkr_result = open('ml_dipole.dat', 'w')
    file_geom = open('geom.xyz', 'r')
    n_atom = int(file_geom.readline())
    file_geom.readline()
    file_kkr_result.write('  ' + str(n_atom) + '\n')
    file_kkr_result.write(' The coordinates\n')
    for i_atom in range(n_atom):
        line_text = file_geom.readline().split()
        line_text[1] = float(line_text[1]) / 0.529
        line_text[2] = float(line_text[2]) / 0.529
        line_text[3] = float(line_text[3]) / 0.529
        file_kkr_result.write(
            str(line_text[0]) + ' ' + str(line_text[1]) + ' ' + str(line_text[2]) + ' ' + str(line_text[3]) + '\n')
    file_geom.close()

    file_kkr_result.write(' Dipole \n')
    for i_y_dim in range(n_y_dim):
        kkr_file = kkr_path + '/' + str(i_y_dim + 1) + '/' + 'gradient.xyz'
        file_kkr_pre = open(kkr_file, 'r')
        file_kkr_pre.readline()
        file_kkr_pre.readline()
        for i_atom in range(n_atom):
            file_kkr_pre.readline()
        line_text = file_kkr_pre.readline().split()
        kkr_energy = energy_zero + float(line_text[3])
        file_kkr_result.write(str(kkr_energy) + '\n')
        file_kkr_pre.close()

    file_kkr_result.write('Gradient of dipole\n')
    for i_y_dim in range(n_y_dim):
        file_kkr_result.write(' State: ' + str(i_y_dim + 1) + '\n')
        kkr_file = kkr_path + '/' + str(i_y_dim + 1) + '/' + 'gradient.xyz'
        file_kkr_pre = open(kkr_file, 'r')
        file_kkr_pre.readline()
        file_kkr_pre.readline()
        for i_atom in range(n_atom):
            file_kkr_pre.readline()

        file_kkr_pre.readline()
        file_kkr_pre.readline()
        file_kkr_pre.readline()

        for i_atom in range(n_atom):
            file_kkr_result.write(file_kkr_pre.readline())
        file_kkr_pre.close()

    file_kkr_result.close()


if __name__ == '__main__':
    start = time.time()
    kkr_single_geom_all()
    elapsed = (time.time() - start)
    print("All time used for each step:",elapsed)
