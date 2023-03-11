# python
import copy
import os
import re
import sys
sys.path.append("../tools/")
import tools


class molpro_log_parser():
    """
    parse molpro log file
    """
    def __init__(self, config={}):
        """ init """
        self.files = {'interface': 'interface.json', 'mo': 'molpro.log'}

        if config != {}:
            root_dir = config['root']
            dirs = config['dirs']
            files = config['files']

            # working directory & files >>>
            self.directory = {}
            self.directory['root'] = root_dir
            self.directory['home'] = root_dir + "/" + dirs['home']
            self.directory[
                'work'] = self.directory['home'] + "/" + dirs['work']

            self.files = {}
            self.files["interface"] = self.directory['work'] + "/" + files[
                'interface']
            self.files[
                "mo"] = self.directory['work'] + "/" + files['molpro_log']
        self.load()
        return

    # -------------------------------------------------------------------------

    def load(self):
        """
        load interface.json
        """
        filename = self.files['interface']
        obj = tools.load_data(filename)
        self.interface = copy.deepcopy(obj)

    def collect_qm(self):
        """
        wrt down in one file
        """
        fileout3 = open('qm_results.dat', 'w')

        qm_interface = tools.load_data("interface.json")
        n_atom = qm_interface['parm']['n_atom']
        geom = qm_interface['mol']
        atoms = geom['atoms']

        fileout3.write('     ' + str(n_atom) + '\n')
        fileout3.write(' The coordinates' + '\n')
        for i in range(n_atom):
            record = atoms[i]
            atomname = record['name']
            coord = record['coord']
            fileout3.write(
                str(atomname) + '   ' + str(coord[0]) + '   ' + str(coord[1]) +
                '   ' + str(coord[2]) + '\n')

        filein4 = open('qm_energy.dat', 'r')
        fileout3.write(filein4.read())
        filein4.close()
# if esa, close
        filein4 = open('qm_gradient.dat', 'r')
        fileout3.write(filein4.read())
        filein4.close()

#        filein4 = open('qm_trdm.dat', 'r')
#        fileout3.write(filein4.read())
#        filein4.close()
# if esa, close
        sourceFile = 'qm_nac.dat'
        if os.path.isfile(sourceFile):
            filein4 = open('qm_nac.dat', 'r')
            fileout3.write(filein4.read())
            filein4.close()
        else:
            for i_state in range(n_state):
                for j_state in range(n_state):
                    fileout3.write('S' + str(i_state) + '    S' +
                                   str(j_state) + '   0.00000   \n')

        sourceFile = 'qm_dipole.dat'
        if os.path.isfile(sourceFile):
            filein4 = open('qm_dipole.dat', 'r')
            fileout3.write(filein4.read())


        fileout3.close()

        return

    # ---------------------------------------------------------------------------

    #   %%% Read the energy
    #   qm_energy.dat
    #   Attention:

    def get_energy(self):
        """ read energy and punch out """

        n_state = self.interface['parm']['n_state']

        logfile = self.files['mo']
        file_in = open(logfile, "r")
        file_out = open("qm_energy.dat", "w")

        file_out.write(' Energy of electronic states' + '\n')

        pattern = re.compile("!MCSCF STATE")

        line = "NOT EMPTY LINE"
        for i_state in range(n_state):
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                    record = line.split()
                    if record[-2] == "Energy":
                        energy = float(record[-1])
                        file_out.write('' + str(energy) + '   ' + '\n')

                        break

        file_in.close()
        file_out.close()

        return

    def get_gradient(self):
        """ read gradient and punch out """

        n_state = self.interface['parm']['n_state']
        n_atom = self.interface['parm']['n_atom']

        logfile = self.files['mo']
        file_in = open(logfile, "r")
        file_out = open("qm_gradient.dat", "w")

        file_out.write(' Gradient of electronic states' + '\n')

        pattern = re.compile("SA-MC GRADIENT FOR STATE")
#        print (file_in)
        line = "NOT EMPTY LINE"
        for i_state in range(n_state):
#        for i_state in range(2):
            file_out.write('  State:           ' + str(i_state + 1) + '\n')
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                    break

            line = file_in.readline()
            line = file_in.readline()
            line = file_in.readline()

            for i_atom in range(n_atom):

                line = file_in.readline()
                record = line.split()

                grad_x = float(record[1])
                grad_y = float(record[2])
                grad_z = float(record[3])

                file_out.write('' + str(grad_x) + '   ' + \
                               str(grad_y) + '   ' + str(grad_z) + '  \n')
        file_in.close()
        file_out.close()

        return

    # -------------------------------------------------------------------------
    #   %%% Read the nac
    #   qm_nac.dat
    #   Attention:
    # ---------------------------------------------------------------------------
    def get_nac(self):

        n_state = self.interface['parm']['n_state']
        n_atom = self.interface['parm']['n_atom']
        n_dime = 3

        nac = []
        for i in range(n_state):
            nac.append([])
            for j in range(n_state):
                nac[i].append([])
                for k in range(n_atom):
                    nac[i][j].append([])
                    for l in range(n_dime):
                        nac[i][j][k].append(0.0)
        """ read nac and punch out """
        logfile = self.files['mo']

        file_in = open(logfile, "r")
#        print (file_in)
        pattern = re.compile("SA-MC NACME FOR STATES")

        line = "NOT EMPTY LINE"
        while line != "":
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                    i = int(line.split()[-3].split('.')[0]) - 1
                    j = int(line.split()[-1].split('.')[0]) - 1
                    break

            line = file_in.readline()
            line = file_in.readline()
            line = file_in.readline()

            k = 0
            while line != "":
                line = file_in.readline()
                if line.strip() == "":
                    break

                record = line.split()

                if (i < n_state) and (j < n_state):

                    nac[i][j][k][0] = -float(record[1])
                    nac[i][j][k][1] = -float(record[2])
                    nac[i][j][k][2] = -float(record[3])

                    nac[j][i][k][0] = -nac[i][j][k][0]
                    nac[j][i][k][1] = -nac[i][j][k][1]
                    nac[j][i][k][2] = -nac[i][j][k][2]

                k = k + 1

        file_in.close()

        file_out = open("qm_nac.dat", "w")

        file_out.write('Nonadiabatic couplings' + '\n')
        for i in range(n_state):
            for j in range(n_state):
                file_out.write(' State:           ' + str(i + 1) + '        ' +
                               str(j + 1) + '\n')
                for k in range(n_atom):
                    file_out.write(
                        str(nac[i][j][k][0]) + '   ' + str(nac[i][j][k][1]) +
                        '    ' + str(nac[i][j][k][2]) + '  \n')

        file_out.close()

        return


    # -------------------------------------------------------------------------
    #   %%% Read the transition diople moments
    #   qm_trdm.dat
    #   Attention:
    # ---------------------------------------------------------------------------

    def get_dipole(self):

        logfile = self.files['mo']

        fileout = open('qm_dipole.dat', 'w')

        filein = open(logfile, "r")

        transition_diploe_moment = []

        self.n_state = self.interface['parm']['n_state']
        self.n_atom = self.interface['parm']['n_atom']

        for i_state in range(self.n_state):
        
           transition_diploe_moment.append([])
        
           for j_state in range(self.n_state):
              
              transition_diploe_moment[i_state].append([])
        
              transition_diploe_moment[i_state][j_state].append(0.0)
              transition_diploe_moment[i_state][j_state].append(0.0)
              transition_diploe_moment[i_state][j_state].append(0.0)
              transition_diploe_moment[i_state][j_state].append(0.0)

        pattern_tdm = re.compile(r'State-averaged charge density matrix saved on record')

        pattern_x = re.compile(r'DMX')
        pattern_y = re.compile(r'DMY')
        pattern_z = re.compile(r'DMZ')

        pattern_finished = re.compile(r'Natural orbital dump')

        pattern_MCSCF = re.compile(r'MCSCF')

        line = "NOT EMPTY LINE"

        while line != "":
        
           line = filein.readline()
        
           tdm_pattern = pattern_tdm.search(line)
        
           if tdm_pattern is not None:
               break

        while line != "":
           line = filein.readline()

           finish_pattern = pattern_finished.search(line)
           if finish_pattern is not None:
               break

           MCSCF_pattern = pattern_MCSCF.search(line)

           if MCSCF_pattern is not None:

              x_pattern = pattern_x.search(line)

              if x_pattern is not None:

                 i_index = int(line[23]) - 1
                 j_index = int(line[31]) - 1

                 record = line.split()
                 transition_diploe_moment[i_index][j_index][0] = float(record[3])
                 transition_diploe_moment[j_index][i_index][0] = transition_diploe_moment[i_index][j_index][0]

              y_pattern = pattern_y.search(line)

              if y_pattern is not None:

                 i_index = int(line[23]) - 1
                 j_index = int(line[31]) - 1

                 record = line.split()
                 transition_diploe_moment[i_index][j_index][1] = float(record[3])
                 transition_diploe_moment[j_index][i_index][1] = transition_diploe_moment[i_index][j_index][1]

              z_pattern = pattern_z.search(line)

              if z_pattern is not None:

                 i_index = int(line[23]) - 1
                 j_index = int(line[31]) - 1

                 record = line.split()
                 transition_diploe_moment[i_index][j_index][2] = float(record[3])
                 transition_diploe_moment[j_index][i_index][2] = transition_diploe_moment[i_index][j_index][2]

        fileout.write('transition dipole moment (X Y Z) in au' + '\n')
        for i in range(self.n_state):
            for j in range(self.n_state):
                fileout.write(' State:           ' + str(i + 1) + '        ' +
                               str(j + 1) + '\n')
                fileout.write(
                        str(transition_diploe_moment[i][j][0]) + '   ' + str(transition_diploe_moment[i][j][1]) +
                        '    ' + str(transition_diploe_moment[i][j][2]) + '  \n')


        filein.close()
        fileout.close()

    def get_other(self):
        """ read transition diople moments and punch out """
#        print("sss")
        n_state = self.interface['parm']['n_state']
        n_atom = self.interface['parm']['n_atom']

        logfile = self.files['mo']
        file_in = open(logfile, "r")
#        file_out = open("qm_trdm.dat", "a")
        file_out = open("qm_other.dat","w")
        file_out.write(' Transition diople moments of electronic states' + '\n')

        pattern = re.compile(r'(!MCSCF trans.+)')

        line = "NOT EMPTY LINE"
#        print (file_in)
        line = file_in.read().splitlines()
        line = [x.strip(' ') for x in line]

        key = []
        for x in line:
            m = pattern.match(x)
            if m:
#       print ('ss')
                line_num = line.index(m.group(1))
                key.append(line[line_num])
        key = key[0:len(key)/2]
        key_x = key[0:len(key)/3]
        key_y = key[len(key)/3:len(key)*2/3]
        key_z = key[len(key)*2/3:len(key)]

        trdm_x = []
        trdm_y = []
        trdm_z = []


        for i in range(n_state):
            for y in key_x:
                n = pattern.match(y)
                if n:
                   trdm_x.append(y)

            for y in key_y:
                n = pattern.match(y)
                if n:
                   trdm_y.append(y)

            for y in key_z:
                n = pattern.match(y)
                if n:
                   trdm_z.append(y)

#        print(trdm_x)
#        print(trdm_y)
#        print(trdm_z)
#       print(trdm_x)
#        for x in range(n_state):
#            for y in range(x+1,n_state):
#                f = x + y*(y - 1)/2
#                file_out.write('{0}   {1}  {2:>18} {3:>18} {4:>18} \n'.format(x,y,trdm_x[f].split()[3],trdm_y[f].split()[3],trdm_z[f].split()[3]))

        file_in.close()
        file_out.close()

        return


    # -------------------------------------------------------------------------
    #   %%% Read all other important information of QM output
    #   molpro.log file is required.
    #   For example: Transition dipole moment and so on
    # ---------------------------------------------------------------------------
    def get_other(self):
        """
        Write other important information in QM output 
        """
        es = []
        gs = []
        pat1e = re.compile("Excited states from <AA,BB:AA,BB> singles matrix")
        pat2e = re.compile("Excitation energies and oscillator strengths")
        float_number = '[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'
        pat1g = re.compile("Charge=(\s)+" + float_number + "(\s)+electrons")
        pat2g = re.compile("XXYZ=(.*)YYXZ=(.*)ZZXY=(.*)")

        # read all
        file_energy = self.files['mo']
        filein = open(file_energy, 'r')

        line = "empty"
        # Excited states from <AA,BB:AA,BB> singles matrix
        while line != "":
            line = filein.readline()
            m1 = pat1e.search(line)
            if m1 is not None:
                break
        line = filein.readline()
        line = filein.readline()

        while line != "":
            line = filein.readline()
            m2 = pat2e.search(line)
            if m2 is not None:
                break
            es.append(line)

            # ground state.
        while line != "":
            line = filein.readline()
            m1 = pat1g.search(line)
            if m1 is not None:
                break
        gs.append(line)
        while line != "":
            line = filein.readline()
            gs.append(line)
            m2 = pat2g.search(line)
            if m2 is not None:
                break
        filein.close()

        fileout = open('qm_other.dat', 'w')
        for line in gs:
            fileout.write(line)
        fileout.write(
            '------------------------------------------------------------- \n')
        for line in es:
            fileout.write(line)
        fileout.write(
            '------------------------------------------------------------- \n')
        fileout.close()

        return


### main program
if __name__ == "__main__":
    ao = molpro_log_parser()

    ao.get_gradient()
    ao.get_nac()
# ao.get_other()
