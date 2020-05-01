#!/usr/bin/env python

# This script was written by manage_regression_suite.py as part of a
# regression_suite file

import sys
import os
import subprocess
import numpy as np
import multiprocessing as mp

os.environ['PYTHONUNBUFFERED'] = '1'
test_failed = False

if not os.path.exists('case_outputs'):
    os.makedirs('case_outputs')

base_path = '/lustre/scratch4/turquoise/.mdt2/jamilg/MPAS-Model/testing_and_setup/compass'


os.chdir(base_path)

locations = []
locations.append('ocean/ha_test/5km/default')
locations.append("ocean/ha_test/10km/default")

commands = []
commands.append(['time', '-p', '/lustre/scratch4/turquoise/.mdt2/jamilg/MPAS-Model/testing_and_setup/compass/ocean/ha_test/5km/default/run_test.py'])
commands.append(['time', '-p', '/lustre/scratch4/turquoise/.mdt2/jamilg/MPAS-Model/testing_and_setup/compass/ocean/ha_test/10km/default/run_test.py'])

datas = []
datas.append([locations[0], commands[0]])
datas.append([locations[1], commands[1]])


def myProcess(data):
    case_output = open('case_outputs/'+data[0].replace("/", "_"), "w")
    print(" ** Running case {}".format(data[0]))
    os.chdir(data[0])
    try:
        subprocess.check_call(data[1], stdout=case_output, stderr=case_output)
        print('      PASS')
    except subprocess.CalledProcessError:
        print('   ** FAIL (See case_outputs/Horizontal_Advection_5km_-_Mesh_Test for more information)')
        test_failed = True
    case_output.close()
    os.chdir(base_path)


p = mp.Pool()
p.map(myProcess, datas)
quit()

print('TEST RUNTIMES:')
case_output = '/case_outputs/'
totaltime = 0
for _, _, files in os.walk(base_path + case_output):
    for afile in sorted(files):
        outputfile = base_path + case_output + afile
        runtime = np.ceil(float(subprocess. check_output(
  ['grep', 'real', outputfile]).decode('utf-8').split('\n')[-2].split()[1]))
        totaltime += runtime
        mins = int(np.floor(runtime/60.0))
        secs = int(np.ceil(runtime - mins*60))
        print('{:02d}:{:02d} {}'.format(mins, secs, afile))
mins = int(np.floor(totaltime/60.0))
secs = int(np.ceil(totaltime - mins*60))
print('Total runtime {:02d}:{:02d}'.format(mins, secs))

if test_failed:
    sys.exit(1)
else:
    sys.exit(0)
