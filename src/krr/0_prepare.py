#!/usr/bin/env python


import os
import sys

import sub_inp_json

sys.dont_write_bytecode = True

filename = 'fitting.input'
xxx = sub_inp_json.read_dat_with_label(filename)
sub_inp_json.dump_json('input_initial.json', xxx)

os.system('cp input_initial.json  ../input_initial.json')
