from glob import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread
import numpy as np


def get_RMSE(folder_name):
  """
  calculates the RMSE of a given resolution
  """

  KPP = glob("../"+folder_name+"/default/forward/output/KPP*")
  KPP = KPP[0]
  print("Processing... {}".format(KPP))


  data = xr.open_dataset(KPP)


  # the center and radius of the cylindars starting point
  simulated_center_location_x =  250000.58 
  radius = 50000
  vi = 1
  nt = 30

  # grab the time in hours, and days
  t1 = str(data.xtime[nt].values)
  sp=t1.find('_')+1
  hr = float(t1[sp:sp+2])
  day = float(t1[sp-3:sp-1])

  # calculate runtime
  totalTime = ((day-1)*24. + hr)*3600.

  # update center
  xLen = data.x_period
  distance = vi * totalTime
  centNew = 250000 + distance

  while centNew > xLen:
    centNew -= xLen

  tracerE = np.zeros_like(data.tracer1[-1,:,:].values)
  last_frame = data.tracer1.shape[0]-1

  for i in range(len(data.xCell.values)):
    if abs(data.xCell.values[i] - centNew) <= radius:
      # if a point is inside the circle set its simulated value to 1
      tracerE[i,:] = 1.0


  # get the RMSE of the simulated values to the actual values
  rmse = np.sqrt(np.mean(tracerE - data.tracer1[-1,:,:].values)**2)

  print("For  {} RMSE = {}".format(folder_name, rmse))

  # return the int(resolution) and the float(rmse)
  return int(folder_name[0:folder_name.find("km")]) , float(rmse)
  

def main():
  folders = ["5km", "10km", "25km"]
  resolution = []
  rmse = []
  for folder in folders:
    x,y = get_RMSE(folder)
    resolution.append(x)
    rmse.append(y)

  print(resolution)
  print(rmse)


  plt.yscale("log")
  plt.xlim(0,26)
  plt.ylim(.0045, .025)
  plt.scatter(resolution,rmse)
  plt.savefig("../visualization/rmse.png")


main()
