#!/usr/bin/env python

import time
import sys
import os
import subprocess
import numpy as np
from multiprocessing import Process
out, err = subprocess.Popen(['nproc'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()

max_procs = int(out.split()[0])
number_of_procs = max_procs

print('Number of usable Processors: {}'.format(number_of_procs))
os.environ['PYTHONUNBUFFERED'] = '1'
test_failed=False
if not os.path.exists('case_outputs'):
    os.makedirs('case_outputs')
base_path = '/lustre/scratch4/turquoise/.mdt2/jamilg/MPAS-Model/testing_and_setup/compass'
os.chdir(base_path)
locations = []
procs = []
commands = []
datas = []


locations.append('ocean/ha_test/5km/default')
commands.append(['time' , '-p' , '/lustre/scratch4/turquoise/.mdt2/jamilg/MPAS-Model/testing_and_setup/compass/ocean/ha_test/5km/default/run_test.py'])
procs.append(8)
datas.append([locations[0], commands[0]])


locations.append('ocean/ha_test/10km/default')
commands.append(['time' , '-p' , '/lustre/scratch4/turquoise/.mdt2/jamilg/MPAS-Model/testing_and_setup/compass/ocean/ha_test/10km/default/run_test.py'])
procs.append(8)
datas.append([locations[1], commands[1]])


locations.append('ocean/ha_test/25km/default')
commands.append(['time' , '-p' , '/lustre/scratch4/turquoise/.mdt2/jamilg/MPAS-Model/testing_and_setup/compass/ocean/ha_test/25km/default/run_test.py'])
procs.append(8)
datas.append([locations[2], commands[2]])


def task(data):
    case_output = open('case_outputs/'+data[0].replace('/', '_'), 'w')
    print(' Running case @ {}'.format(data[0]))
    print(' Running command: {}'.format(data[1]))
    os.chdir(data[0])
    subprocess.check_call(data[1])

start_time = time.time()

running = []
index = 0

while True:

  assert max_procs >= number_of_procs

  print("current index : {}".format(index))
  print("number_of_procs processes: {}".format(number_of_procs))
  if index > len(procs)-1:
    break ;

  if number_of_procs >= procs[index]:
    print("we have {} number_of_procs, current task uses {} added to queue".format(number_of_procs, procs[index]))
    number_of_procs = number_of_procs - procs[index]
    proc = Process(target=task, args=(datas[index],))
    running.append([proc, procs[index]])
    proc.start()
    
  
  if index == len(procs)-1 or  number_of_procs < procs[index] :
    print("not enough procs waiting till some finish or started all in queue")
    while True:
      if number_of_procs >= procs[index]:
        break
      else:
        for proc in running:
          proc[0].join()
          if not proc[0].is_alive():
            number_of_procs += proc[1]
            running.remove(proc)

   
  index = index + 1   



end_time = time.time()
print('parallel run time: {} min'.format((end_time - start_time) / 60))
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
