#!/usr/bin/env python

from math import *
import sys
import numpy
import raytrace

def vmag(vect):
    return sqrt(numpy.dot(vect,vect))

def vnorm(vect):
    return vect/vmag(vect)

narg = len(sys.argv)

x = []
y = []
z = []
kx = []
ky = []
kz = []
n = []
z_vert = [0.0]
R = [0.0]
K = [0.0]
in_ray = [0.0]*3
out_ray = [0.0]*3
surf_params = [0.0]*5

screen_pos=0.0

if (narg > 7):			# must have at least these four
    filename = sys.argv[1]
    x.append(float(sys.argv[2]))
    y.append(float(sys.argv[3]))
    z.append(float(sys.argv[4]))
    kx.append(float(sys.argv[5]))
    ky.append(float(sys.argv[6]))
    kz.append(float(sys.argv[7]))
else:
    print "Must supply lens_file_name x0, y0, z0, kx0, ky0, kz0 arguments"
    sys.exit()

if (narg > 8):		# optionally, put a screen somewhere
    screen_pos = float(sys.argv[8])

lens_file = open(filename,'r');	# grab lens surface parameters
n_surf = int(lens_file.readline().strip())	# number of surfaces (1st line)
n.append(float(lens_file.readline().strip()))	# initial refr. index (2nd line)
current_z = 0.0;
for i in range(n_surf): 		# and n_surf additional lines...
    # read in lens file and verify results
    line = lens_file.readline()
    Slist = line.split()
    n.append(float(Slist[0]))
    z_vert.append(float(Slist[1]))
    R.append(float(Slist[2]))
    K.append(float(Slist[3]))
    current_z += z_vert[i+1]
    print "Surface %d has n = %f, z_vert = %f, radius = %g, K = %f" % \
                 (i+1,n[i+1],current_z,R[i+1],K[i+1])

lens_file.close()

print "Ray 1 has x, k = (%f,%f,%f), (%f,%f,%f)" % \
                (x[0],y[0],z[0],kx[0],ky[0],kz[0])

current_z = 0.0
for i in range(n_surf):		# now propagate surface-at-a-time
				# begin ray propagation for loop
    # populate surface parameters array
    current_z += z_vert[i+1]
    surf_params[0] = n[i]
    surf_params[1] = current_z
    surf_params[2] = R[i+1]
    surf_params[3] = K[i+1]
    surf_params[4] = n[i+1]

    # populate in_ray array for +y ray
    ray_pos = numpy.array([x[i],y[i],z[i]],dtype='d')
    ray_direc = vnorm(numpy.array([kx[i],ky[i],kz[i]],dtype='d'))

    # carry out calculation
    out_pos,out_direc,back = raytrace.ray_step_3d(ray_pos,ray_direc,surf_params)
    if back:
        print "WARNING: Ray jumped backwards"

    # stow out_ray into approp. arrays
    x.append(out_pos[0])
    y.append(out_pos[1])
    z.append(out_pos[2])
    kx.append(out_direc[0])
    ky.append(out_direc[1])
    kz.append(out_direc[2])

    print "Ray %d has x, k = (%f,%f,%f), (%f,%f,%f)" % \
                (i+2,x[i+1],y[i+1],z[i+1],kx[i+1],ky[i+1],kz[i+1])

# compute intercepts
s_screen = (screen_pos - z[n_surf])/kz[n_surf]
x_screen = x[n_surf] + s_screen*kx[n_surf]
y_screen = y[n_surf] + s_screen*ky[n_surf]

s_zy = -x[n_surf]/kx[n_surf]
y_zy = y[n_surf] + s_zy*ky[n_surf]
z_zy = z[n_surf] + s_zy*kz[n_surf]

s_zx = -y[n_surf]/ky[n_surf]
x_zx = x[n_surf] + s_zx*kx[n_surf]
z_zx = z[n_surf] + s_zx*kz[n_surf]

print "Ray intercepts screen at (%f, %f, %f)" % (x_screen,y_screen,screen_pos)
print "Ray intercepts z-y plane at (%f, %f, %f)" % (0.0,y_zy,z_zy)
print "Ray intercepts z-x plane at (%f, %f, %f)" % (x_zx,0.0,z_zx)

