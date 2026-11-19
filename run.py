#!/usr/bin/python3

# SPDX-FileCopyrightText: Copyright 2023, 2025 Harusato Kimura
# SPDX-License-Identifier: MIT

# ==============================================================================

import sys
import os
import shutil
import json
import math
import subprocess
import numpy as np
import argparse

# ==============================================================================

parser = argparse.ArgumentParser()

# I/O
parser.add_argument("-if", "--input-file", type=str, required=True)
parser.add_argument("-od", "--output-directory", type=str)
# Calculation params
parser.add_argument("-nfp", "--num-frequency-points", type=int, default=250)
parser.add_argument("-ubf", "--upper-bound_frequency", type=float, default=25.0)
parser.add_argument("-ngp", "--num-grid-points", type=int, default=100)

args = parser.parse_args()

nf = args.num_frequency_points
fmax = args.upper_bound_frequency
ngrid = args.num_grid_points
df = fmax / float(nf + 1)


path2paramfile = args.input_file
targ_dir = os.path.dirname(path2paramfile) + "/"
case_name = os.path.splitext(os.path.basename(path2paramfile))[0]

res_dir = targ_dir + case_name + "_results_grtm/"
if args.output_directory is not None:
    res_dir = args.output_directory
if res_dir[-1] != "/":
    res_dir = res_dir + "/"

if os.path.exists(res_dir):
    shutil.rmtree(res_dir)
os.makedirs(res_dir)

# ==============================================================================

f = open(res_dir + "/employed_params.txt", "w")
f.write("nf = {} points\n".format(nf))
f.write("fmax = {} Hz\n".format(fmax))
f.write("df = {} Hz\n".format(df))
f.write("ngrid = {} points\n".format(ngrid))
f.close()

# ==============================================================================

argv = " " + path2paramfile
argv += " " + res_dir
argv += " " + str(nf)
argv += " " + str(fmax)
argv += " " + str(ngrid)
subprocess.run(["./bin/main" + argv], shell=True)
