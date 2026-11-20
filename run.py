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
parser.add_argument(
    "-nfp",
    "--num-frequency-points",
    type=int,
    default=250,
    help="Number of frequency intervals nf. The frequencies are "
    "fmin + (fmax - fmin) / nf * i for i = 0 to nf (nf + 1 points).",
)
parser.add_argument(
    "-lbf",
    "--lower-bound-frequency",
    type=float,
    default=0.0,
    help="Lower bound of frequency fmin (Hz).",
)
parser.add_argument(
    "-ubf",
    "--upper-bound-frequency",
    "--upper-bound_frequency",
    dest="upper_bound_frequency",
    type=float,
    default=25.0,
    help="Upper bound of frequency fmax (Hz).",
)
parser.add_argument(
    "-ngp",
    "--num-grid-points",
    type=int,
    default=100,
    help="Number of grid points per predicted mode for the root search "
    "along the phase velocity.",
)
# Only for a profile with Qp and Qs (viscoelastic medium)
parser.add_argument(
    "-el",
    "--elastic",
    action="store_true",
    help="Ignore Qp and Qs of the profile and treat the medium as elastic "
    "(the medium responses and the relative powers are computed instead "
    "of the attenuation coefficients).",
)
parser.add_argument(
    "-f0",
    "--reference-frequency",
    type=float,
    default=None,
    help="Reference frequency (Hz) at which Vp and Vs of the profile are "
    "defined. If given, the causal constant-Q model of Kjartansson (1979) "
    "is used; otherwise the complex velocities are frequency independent.",
)

args = parser.parse_args()

nf = args.num_frequency_points
fmin = args.lower_bound_frequency
fmax = args.upper_bound_frequency
ngrid = args.num_grid_points
f_ref = args.reference_frequency

if nf < 2:
    parser.error("--num-frequency-points must be 2 or larger")
if fmin < 0.0:
    parser.error("--lower-bound-frequency must not be negative")
if fmax <= fmin:
    parser.error("--upper-bound-frequency must be larger than --lower-bound-frequency")
if ngrid < 1:
    parser.error("--num-grid-points must be positive")
if f_ref is not None and f_ref <= 0.0:
    parser.error("--reference-frequency must be positive")

df = (fmax - fmin) / float(nf)


path2paramfile = args.input_file
targ_dir = os.path.dirname(path2paramfile) + "/"
case_name = os.path.splitext(os.path.basename(path2paramfile))[0]

res_dir = targ_dir + case_name + "_results_grtm/"
if args.output_directory is not None:
    res_dir = args.output_directory
if res_dir[-1] != "/":
    res_dir = res_dir + "/"

# A profile with Qp and Qs (6 values per layer) is a viscoelastic medium, for
# which the attenuation coefficients are computed instead of the medium
# responses and the relative powers.
with open(path2paramfile) as f:
    layers = [l for l in f if l.strip() != "" and not l.strip().startswith("#")]
attenuation_given = len(layers) > 0 and len(layers[0].split(",")) == 6
viscoelastic = attenuation_given and not args.elastic
if f_ref is not None and args.elastic:
    parser.error("--reference-frequency cannot be used with --elastic")
if f_ref is not None and not viscoelastic:
    parser.error("--reference-frequency requires a profile with Qp and Qs")

if os.path.exists(res_dir):
    shutil.rmtree(res_dir)
os.makedirs(res_dir)

# ==============================================================================

f = open(res_dir + "/employed_params.txt", "w")
# nf is the number of frequency intervals (nf + 1 points from fmin to fmax).
f.write("nf = {} points\n".format(nf))
f.write("fmin = {} Hz\n".format(fmin))
f.write("fmax = {} Hz\n".format(fmax))
f.write("df = {} Hz\n".format(df))
f.write("ngrid = {} points\n".format(ngrid))
if attenuation_given and not viscoelastic:
    f.write("medium = elastic (Qp and Qs of the profile were ignored)\n")
if viscoelastic:
    if f_ref is None:
        f.write("velocity_model = frequency independent\n")
        f.write("reference_frequency = None\n")
    else:
        f.write("velocity_model = Kjartansson (1979)\n")
        f.write("reference_frequency = {} Hz\n".format(f_ref))
f.close()

# ==============================================================================

argv = " " + path2paramfile
argv += " " + res_dir
argv += " " + str(nf)
argv += " " + str(fmax)
argv += " " + str(ngrid)
argv += " " + str(fmin)
argv += " " + ("none" if f_ref is None else str(f_ref))
argv += " " + ("elastic" if args.elastic else "auto")
result = subprocess.run(["./bin/main" + argv], shell=True)
sys.exit(result.returncode)
