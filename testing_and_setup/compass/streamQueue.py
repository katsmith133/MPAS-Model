from multiprocessing import Process
import time
import subprocess
available = 40
tasks = [10,10,20,20,10]

running = []
index = 0

def task(message):
  time.sleep(1)
  single_proc = subprocess.Popen([message],shell=True,stdout=subprocess.PIPE)
  print(single_proc.communicate()[0])

while True:
  print("current index : {}".format(index))
  print("available processes: {}".format(available))
  if index >= len(tasks):
    break ;

  if available >= tasks[index]:
    print("we have {} available, current task uses {} added to queue".format(available, tasks[index]))
    available = available - tasks[index]
    command = "echo available {}, used {}".format(available, tasks[index])
    proc = Process(target=task, args=(command,))
    running.append([proc, tasks[index]])
    proc.start()
    index = index + 1
  else:
    print("not enough procs waiting till some finish")
    print("\tavailable: {}\n\tindex: {}\n\tnext task: {}".format(available, index, tasks[index]))
    while True:
      if available >= tasks[index]:
        break
      else:
        for proc in running:
          proc[0].join(timeout=0)
          if not proc[0].is_alive():
            available += proc[1]

