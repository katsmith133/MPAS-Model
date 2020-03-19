#!/usr/bin/env python


# Read in the data
import xarray as xr
# get files
import glob
# get unique
import numpy as np
# plot
import matplotlib.pyplot as plt
def main():
  filename = None
  filename = glob.glob("forward/output/KPP*")

  if filename == None:
    print("No output files found")
    quit()
  else: 
    filename = filename[0]

  data = xr.open_dataset(filename)

  xCell_unique = np.unique(data.xCell.values)
  zMid_unique  = np.unique(data.zMid.values)
  tracer1 = data.tracer1.values[0,:100,:]



  plt.contourf(xCell_unique[:100], zMid_unique[:100], tracer1)
  plt.savefig("t.png")

main()
