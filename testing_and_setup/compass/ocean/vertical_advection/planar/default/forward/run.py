#!/usr/bin/env python

# This script was generated from setup_testcases.py as part of a config file

import sys
import os
import shutil
import glob
import subprocess


dev_null = open('/dev/null', 'w')

# Run command is:
# gpmetis graph.info 8
subprocess.check_call(['gpmetis', 'graph.info', '8'])

# Run command is:
# mpirun -n 8 ./ocean_model -n namelist.ocean -s streams.ocean
subprocess.check_call(['mpirun', '-n', '8', './ocean_model', '-n',
                       'namelist.ocean', '-s', 'streams.ocean'])
