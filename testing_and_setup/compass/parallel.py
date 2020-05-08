#!/usr/bin/env python

import time
import sys
import os
import subprocess
import psutil
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

start_time = time.time()

Queue_running = []
index = 0
continue_add = True
Done = False
print_index = 0
base = os.getcwd()
while True:
  if continue_add and not Done:
    print("-----------IN ADDING PHASE----------")
    if number_of_procs >= procs[index]:
      print("we have {} number_of_procs, current task uses {}\n\tAdding to queue".format(number_of_procs, procs[index]))
      try:
        case_output = open('case_outputs/'+datas[index][0].replace('/', '_'), 'w+')
        #chdir_proc = "cd " + str(base) + "/"+str(datas[index][0]) + ";"
        #chdir_base = "cd " + str(base) + ";"
        #datas[index][1].insert(0, chdir_proc)
        #datas[index][1].append(chdir_base)
        print("processing command: {}".format(datas[index][1])) 
        os.chdir(datas[index][0])
        open_proc = subprocess.Popen(datas[index][1],  stdout=case_output, stderr=case_output)
        os.chdir(base)
        Queue_running.append([open_proc, procs[index], datas[index][1]])
        number_of_procs = number_of_procs - procs[index]
        print("\tNew number_of_procs: {}".format(number_of_procs))
      except subprocess.CalledProcessError:
        print("error for : " + str(data[1]))
    elif number_of_procs < procs[index]:
      print("NOT ENOUGH WAIT FOR PROCESSES")
      continue_add = False
    index = index + 1
    if index > len(procs) -1:
      print("No more to add moving to processing phase")
      continue_add = False

  elif not continue_add and not Done:
    if print_index % 10000000 == 0:
      print("Cant add and not done")
    for background_process in Queue_running:
      print_index = print_index + 1
      background_process[0].wait()
      pid = background_process[0].pid
      if not psutil.pid_exists(pid):
        Queue_running.remove(background_process)
        print(str(pid) +" compleated")
        number_of_procs = number_of_procs + background_process[1]
        if index < len(procs)-1:
          print("STILL HAVE MORE TO PROCESS")
          if number_of_procs >= procs[index]:
            print("enough processor free; going to add stage")
            continue_add = True
            #break
        else:
          print("NO MORE TO ADD")
          if len(Queue_running) == 0:
            print("all running processes Done leaving")
            Done = True
            continue_add = False
  elif Done and not continue_add:
    print("done with all")
    break 


end_time = time.time()
print('parallel run time: {} min'.format((end_time - start_time) / 60))

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
