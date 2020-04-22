Vertical advection test
-------------------------

This vertical advection test has the following parameters:


Parameter definition
--------------------
  * Planar hex:
    * dc, nx , ny: values for generating a mesh 
  * Sphere parameters:
    * Radius of the sphere
    * x_cent , y_cent
    * vertVelocityTopVar: ??speed moving in the Z direction??
    * ZonalC: speed moving in the x direction
    * MeridionalC: speed moving in the y direction
  * layers:
    * layer_1: ???
    * layer_2: ???
  * run_duration: simulation run duration


Parameter Values
-------------------
  * Planar hex:
    * dc, nx , ny = 100
  * Sphere parameters:
    * Radius of the sphere = 100
    * x_cent , y_cent = 5000, 5000
    * vertVelocityTopVar = .05 
    * ZonalC = 0
    * MeridionalC = 0
  * layers:
    * layer_1 = 70
    * layer_2 = 90


Explanation
----------------------
The grid size is 10000 x 8660.25403; this can be calcualted as follows
  * x = dc * nx
  * y = [ sqrt(3)/2 ] * (dc * ny)

Thus:
  * x = 100 * 100 = 10000
  * y = [ sqrt(3) / 2 ] * (100 * 100) = 8660.25


