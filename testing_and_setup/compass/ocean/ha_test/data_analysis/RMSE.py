from glob import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread
import numpy as np

def is_inside(point, center, radius):
  """

  checks a given point to see if it is inside
  the circle at center with radius

  point  = 2D list with x,y coordinant
  center = 2D list with x,y coordinants for center of circle
  radius = the lenght of the radius for the circle

  """
  # grab X and Y values for point and center of circle
  x = point[0]
  y = point[1]

  center_x = center[0]
  center_y = center[1]


  # calculates the distance from the point to the center
  is_inside = abs(x - center_x)

  # if the distance is <= radius; it is inside the circle
  return is_inside <= float(radius)



def get_RMSE(folder_name):
  """
  calculates the RMSE of a given resolution
  """
  # the center and radius of the cylindars starting point
  simulated_center_location_x =  250000.58 
  simulated_center_location_y =  216506.35 
  radius = 50000

  # cylanders inital velocity
  vi = 1

  # The offset (in seconds) from running an int number of days
  # seconds = 86400 * (runtime days - actual time days )
  offset_time = 86400 * ( 6 - 5.787037037037037 )

  # calculating the distance traveled
  distance = vi * offset_time

  # new x position
  simulated_center_location_x = simulated_center_location_x + distance

  # generate the circle

  # Grab the KPP file from the specific resolution study
  KPP = glob("../"+folder_name+"/default/forward/output/KPP*")
  KPP = KPP[0]
  print("Processing... {}".format(KPP))


  # load data
  data = xr.open_dataset(KPP)
  last_frame = data.tracer1.shape[0]-1

  print(last_frame)
  # grab the number of cells
  nCells = data.dims['nCells']
  tracer_exact = np.zeros(nCells)

  
  # get points inside location
  index = 0
  inside_x = []
  inside_y = []

  for i in zip(data.xCell.values, data.yCell.values):
    if is_inside(i, [simulated_center_location_x,simulated_center_location_y], radius):
      # if a point is inside the circle set its simulated value to 1
      inside_x.append(i[0])
      inside_y.append(i[1])
      tracer_exact[index] = 1
  index = index + 1

  # get the RMSE of the simulated values to the actual values
  rmse = np.sqrt(np.mean(tracer_exact - data.tracer1[last_frame,:,99].values)**2)

  print("For  {} RMSE = {}".format(folder_name, rmse))

  # plot the simulated circle over the actual values
  plt.scatter(data.xCell, data.yCell, c=data.tracer1[last_frame,:,99], vmin=0, vmax=1)
  plt.scatter(inside_x, inside_y,c="RED", alpha=.3)
  plt.savefig("../visualization/"+str(folder_name)+"_plot.png")

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
  plt.ylim(.17, .2)
  plt.scatter(resolution,rmse)
  plt.savefig("../visualization/rmse.png")


main()
