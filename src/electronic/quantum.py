#! /usr/bin/env python2

from tools.namelist import *
from interface_molpro_qm.molpro import *

import os
import sys
import shutil
import time

sys.dont_write_bytecode = True


# job wrapper of different quantum-chemistry package.

# % read fortran namelist to get md info. & the target qc package. 
#   convert to json format.
# % call related quantum code
#
# %
#

class QuanChem:
    """
    model for quantum chemistry calculations.
    """

    def __init__(self):

        self.files = {}
        self.files["interface"] = "qm_interface"
        self.files["dyn"] = "dyn.inp"
        self.files['dyn_json'] = "inp.json"

        return

    def prepare(self):
        """
        read md info. & prepare for the qc job
        """
        # interface file
        # generate interface.json
        ic = interface_converter(filename=self.files['interface'])
        i_time = int(ic['parm']['i_time'])

        if i_time == 0:
            # read dynamic in file
            print "zero step"
            nma = namelist(filename=self.files['dyn'])
            self.obj = nma.get()
        else:
            nma = namelist(filename=self.files['dyn'])
            self.obj = nma.get()
#            print "STEP : ", i_time

        return

    def run(self):
        """ 
        distributing the job
        """
        qm_package = int(self.obj['quantum']['qm_package'])

        if qm_package == 105:
            # call molpro
#            print "molpro is running"
            Molpro()
        else:
            # not done.
            # exit.
            print "cannot work with this qm package type: ", qm_package
            sys.exit(1)

        return

    def finalize(self):
        """ action after qc. """

        print "ELECTRONIC CALCULATION JOB DONE"

        return


if __name__ == "__main__":
    start_all = time.time()
    qc = QuanChem()
    qc.prepare()
    qc.run()
    qc.finalize()
    elapsed_all = (time.time() - start_all)
    print("Time used for electronical structure calculation",elapsed_all)
