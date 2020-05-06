from multiprocessing import Process
import time
import subprocess
max_procs = 40
number_of_procs = 40
tasks = [10,10,20,20,10]

running = []
index = 0

def task(message):
  time.sleep(1)
  single_proc = subprocess.Popen([message],shell=True,stdout=subprocess.PIPE)
  print(single_proc.communicate()[0])

while True:
  assert max_procs >= number_of_procs
  print("current index : {}".format(index))
  print("number_of_procs processes: {}".format(number_of_procs))
  if index >= len(tasks):
    break ;

  if number_of_procs >= tasks[index]:
    print("we have {} number_of_procs, current task uses {} added to queue".format(number_of_procs, tasks[index]))
    number_of_procs = number_of_procs - tasks[index]
    command = "echo number_of_procs {}, used {}".format(number_of_procs, tasks[index])
    proc = Process(target=task, args=(command,))
    running.append([proc, tasks[index]])
    proc.start()
    index = index + 1
  else:
    print("not enough procs waiting till some finish")
    print("\tnumber_of_procs: {}\n\tindex: {}\n\tnext task: {}".format(number_of_procs, index, tasks[index]))
    while True:
      if number_of_procs >= tasks[index]:
        break
      else:
        for proc in running:
          proc[0].join(timeout=0)
          if not proc[0].is_alive():
            number_of_procs += proc[1]

