
Nx = int(input("enter Nx:  "))
Ny = int(input("enter Ny:  "))
dc = int(input("enter dc:  ")) 
vi = float(input("enter init velocity:  "))

# length of x
x = dc * Nx
# length of y
y = ( (3**.5)/2 ) * dc * Ny

# we will start in the middle, taking up 20% of the screen
Init_x_point = round( x/2 , 2 )
Init_y_point = round( y/2 , 2 )
radius = round(.1 * x, 2)

# t = d / vi 
runtime = round(x  / 86400 , 2 ) 
# distance traveled = d = t * vi
# displacement = inital_point + distrance traveled
dis_travel = round(Init_x_point + (runtime * vi) , 2)

# back in domain = distance traveled % domain size
displace = round(dis_travel % x , 2)

print("\n\n\n\nNx = {}\n Ny = {}\n dc = {}\n vi = {}\n x lenght = {}\n y length = {}".format(Nx, Ny, dc, vi, x, y))
print("\ninit_point = ( {} , {} )".format(Init_x_point, Init_y_point))
print("\nradius = {}\n runtime = {} days".format(radius, runtime ))

print("new location = ( {}, {} )".format(displace, Init_y_point))
