#!/usr/bin/env python

# This script was generated from setup_testcases.py as part of a config file

import sys
import os
import shutil
import glob
import subprocess


dev_null = open('/dev/null', 'w')

# Run command is:
# MpasMeshConverter.x base_mesh.nc mesh.nc
subprocess.check_call(['MpasMeshConverter.x', 'base_mesh.nc', 'mesh.nc'])

# Run command is:
# mpirun -n 1 ./ocean_model -n namelist.ocean -s streams.ocean
subprocess.check_call(['mpirun', '-n', '1', './ocean_model', '-n',
                       'namelist.ocean', '-s', 'streams.ocean'])
