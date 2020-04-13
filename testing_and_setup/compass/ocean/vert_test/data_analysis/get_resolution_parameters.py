#####################################
#
#
# This is for a 25k vertical 
# resolution study
#
# If I wanted this to be for a
# N resolution study
# I would need to ...
#
#
####################################



# -----------variables---------------
# doesnt change for vertical advection
Nx = 20

# doesnt change for vertical advection
Ny = 20

# doesnt change for vertical advection
dc = 25000 

# I want to keep this the same
#????If this is going from lvl 5 to bottom shouldnt it be negitive????
vi = .05

# based off user input (ie, 50, 25, 10, 5..)
dz =  int(input("enter DZ (m): ")) 

# this doesnt change for vertical advection
bottom_depth = 1000

# calculated
vert_levels = bottom_depth / dz



# --------calculations-------------
# length of x
x = dc * Nx
# length of y
y = ( (3**.5)/2 ) * dc * Ny

# we will start in the middle, taking up 10% of the screen
Init_x_point = round( x/2 , 2 )
Init_y_point = round( y/2 , 2 )
Init_z_point = bottom_depth - 50

radius = round(.1 * x, 2)

# t = d / vi 
# we only want to travel 3000m not length of x
runtime = round( (bottom_depth - 100) / vi , 2 )
 
# distance traveled = d = t * vi
dis_travel = round(Init_z_point + (runtime * vi) , 2)

#---------------print results------------


print("--nx: {}\n--ny: {}\n--dc: {}\n".format(Nx, Ny, dc))

print("Location:\n\tx_cent: {}\n\ty_cent: {}\n\t??z??: {}".format(Init_x_point, Init_y_point, Init_z_point))


print("\n??dz??: {}\nconfig_vertical_advection_bottom_depth: {}\n??vert_levels??: {}".format(dz, bottom_depth,vert_levels))

print("layer_1: {}\nlayer_2: {}".format("--?--", "--?--"))

print("\nconfig_tracer_radius: {}".format(radius))

print("config_run_duration: {} days\n".format(runtime/86400 ))

print("distance_travled = {} ".format(dis_travel))



print("\n\n\ntracer will start in the z-axis @ {} and move up at a rate of {}\n its will run for {} days and stop at {}".format(Init_z_point, vi, runtime/86400, dis_travel))
