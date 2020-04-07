Horizontal advection test
-------------------------

This horizontal advection test has the following parameters:

  * Planar hex: dc = nx = ny = 100
  * Radius of sphere = 2000
  * Starting point of sphere = (5000, 7500) 
  * Merininal velocity = 0
  * Zonal velocity = .1


The inital location of the sphere is meant to overflow towards the bottom of the mesh and
move horizontally. It should be shown that as the sphere hits the left most boundry it
will overflow to the right boundry. Since this is the horizontal advection we set the 
Merinial vecloity to zero.

The grid size is 10000 x 8660.25403; this can be calcualted as follows
  * x = dc * nx 
  * y = [ sqrt(3)/2 ] * (dc * ny)

Thus:
  * x = 100 * 100 
  * y = [ sqrt(3) / 2 ] * (100 * 100)

The cylander travels .1 meter per second for 86400 seconds; thus the cylander
travels 8640 meters per day. This is running for 1 day thus it travels 8640 meter

Original center is (5000, 7500) ; the new location is
  * 5000 + 8640 = (13640 , 7500)
  * This is outside the domain - push inside the domain (13640 % 10000) = 3640
  * New location = (3640 , 7500)


with new x, y, rebuild cylinder based on radius
xCell and yCell from MPAS

tracer1_exact(:,:) = 0.0
do i = 1, nCells
   dist = sqrt((xCell[i] - x_new)**2 + (yCell[i] - y_new)**2)
   if dist < radius
      tracer1_exact(:,i) = 1
enddo



do the RMSE with MPAS data
