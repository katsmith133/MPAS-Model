#!/usr/bin/env python


# Read in the data
import xarray as xr
# get files
from glob import glob
# get unique
import numpy as np
# plot
import matplotlib.pyplot as plt


def get_RMSE(folder_name):
  KPP = glob("../"+folder_name+"/default/forward/output/KPP*")
  KPP = KPP[0]
  print("Processing... {}".format(KPP))


  data = xr.open_dataset(KPP)

  xCell_unique = np.unique(data.xCell.values)
  zMid_unique  = np.unique(data.zMid.values)
  tracer1 = data.tracer1.values[0,:100,:]


  
  bottomDepth = data.bottomDepth.values
  config_vertical_advection_layer_thickness = data.config_vertical_advection_layer_thickness

  nt = 0
  # grab the time in hours, and days
  t1 = str(data.xtime[nt].values)
  sp=t1.find('_')+1
  hr = float(t1[sp:sp+2])
  day = float(t1[sp-3:sp-1])

  # calculate runtime
  totalTime = ((day-1)*24. + hr)*3600.
  
  velocity = 1
  distance_traveled = totalTime * velocity

  topofBox = bottomDepth[0] - config_vertical_advection_layer_thickness - distance_traveled
  bottomofBox = bottomDepth[0] - distance_traveled


  tracerE = np.zeros_like(data.tracer1[-1,:,:].values)
  if folder_name == "5m":
    print("\n\nbottomDepth: {}".format(bottomDepth.shape))
    print("topofBox:    {}".format(topofBox))
    print("bottomofBox: {}".format(bottomofBox.shape))
    print("zMidUnique:  {}\n\n\n".format(zMid_unique.shape))


  for k in range(len(zMid_unique)):
    if zMid_unique[k] >= topofBox:
      if zMid_unique[k] <= bottomofBox:
        tracerE[k,:] = 1.0
  
  rmse = np.sqrt(np.mean(tracerE - data.tracer1[-1,:,:].values)**2)
  print("For  {} RMSE = {}".format(folder_name, rmse))


  return int(folder_name[0:folder_name.find("m")]) , float(rmse)


def main():
  folders = ["5m", "10m", "20m"]#,"25m", "50m","100m"]

  resolution = []
  rmse = []
  for folder in folders:
    x,y = get_RMSE(folder)
    resolution.append(x)
    rmse.append(y)

  print(resolution)
  print(rmse)

  # find line of best fit
  resolution = np.array(resolution)
  rmse = np.log10(np.array(rmse))
  m,b = np.polyfit(resolution, rmse,1)

  points = np.linspace(min(resolution), max(resolution),100)
  plt.plot(points, m*points+b)
  plt.title("y = {} x + {}".format(m,b))
  plt.scatter(resolution,rmse)
  plt.savefig("../visualization/rmse.png")


main()
