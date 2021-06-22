#!/usr/bin/env python

from math import *
import sys

def    ray_step(in_ray, surface):

    out_ray = [0.0]*3

    # rename in_ray components for clarity
    z0 = in_ray[0]
    y0 = in_ray[1]
    thet0 = in_ray[2]

    # rename surface params for clarity
    n0 = surface[0]
    z_v = surface[1]
    R = surface[2]
    K = surface[3]
    n_new = surface[4]

    # establish sign variable for sign flips
    Rsign = R/fabs(R)
    Ksign = 1.0
    if (K < -1.0):
        Ksign *= -1.0

    # compute solution to quadratic formula
    m = tan(thet0)
    if (K == -1.0 and m == 0.0):	# would make denom = 0
        z_new = z_v + y0*y0/(2.0*R)	# parabolic case
    else:
        denom = 1.0 + K + m*m
        b = (m*y0 - m*m*z0 - R - (K+1)*z_v)/denom
        c = (y0*y0 + m*m*z0*z0 - 2.0*m*z0*y0 + 2.0*R*z_v + (K+1)*z_v*z_v)/denom
        z_new = -b - Rsign*Ksign*sqrt(b*b - c)
    if (fabs(R) > 1.0e10):
        z_new = z_v	# planar case

    y_new = y0 + m*(z_new - z0)

    thet_norm = atan(-y_new*Rsign/sqrt(R*R - (K+1)*y_new*y_new))
    thet_in = thet0 - thet_norm
    thet_out = asin(n0*sin(thet_in)/n_new)
    thet_new = thet_out + thet_norm
    #print("%f %f %f %f" % (z_v,z_new,y_new,thet_norm))

    if (z_new < z0):
        print("WARNING: ray jumped backwards!")

    out_ray[0] = z_new
    out_ray[1] = y_new
    out_ray[2] = thet_new

    return out_ray

narg = len(sys.argv)

z = []
y = []
thet = []
n = []
z_vert = [0.0]
R = [0.0]
K = [0.0]
in_ray = [0.0]*3
out_ray = [0.0]*3
surf_params = [0.0]*5

screen_pos=0.0

if (narg > 5):			# must have at least these five
    filename = sys.argv[1]
    startz = float(sys.argv[2])
    starty = float(sys.argv[3])
    slope = float(sys.argv[4])
    offset = float(sys.argv[5])
else:
    print("Must supply lens_file_name z0, y0, slope, offset arguments")
    sys.exit()

if (narg > 6):		# optionally, put a screen somewhere
    screen_pos = float(sys.argv[6])

# define +y ray to be one offset above start position
z.append(startz)
y.append(starty + offset)
thet.append(atan(slope))	# input slope --> angle in radians

# define -y ray to be one offset below start position; same angle
z.append(startz)
y.append(starty - offset)
thet.append(thet[0])

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
    print("Surface %d has n = %f, z_vert = %f, radius = %g, K = %f" % \
                 (i+1,n[i+1],current_z,R[i+1],K[i+1]))

lens_file.close()

print("Ray 1+ has z = %f; y = %f; thet = %f" % (z[0],y[0],thet[0]))
print("Ray 1- has z = %f; y = %f; thet = %f" % (z[1],y[1],thet[1]))

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
    in_ray[0] = z[2*i]
    in_ray[1] = y[2*i]
    in_ray[2] = thet[2*i]

    # carry out +y calculation
    out_ray = ray_step(in_ray,surf_params)

    # stow out_ray into approp. arrays
    z.append(out_ray[0])
    y.append(out_ray[1])
    thet.append(out_ray[2])

    # populate in_ray array for -y ray
    in_ray[0] = z[2*i+1]
    in_ray[1] = y[2*i+1]
    in_ray[2] = thet[2*i+1]

    # carry out -y calculation
    out_ray = ray_step(in_ray,surf_params)

    # stow out_ray into approp. arrays
    z.append(out_ray[0])
    y.append(out_ray[1])
    thet.append(out_ray[2])

    print("Ray %d+ has z = %f; y = %f; thet = %f" % \
                            (i+2,z[2*i+2],y[2*i+2],thet[2*i+2]))
    print("Ray %d- has z = %f; y = %f; thet = %f" % \
                            (i+2,z[2*i+3],y[2*i+3],thet[2*i+3]))

# final +y ray parameters
zf0 = z[2*n_surf]
yf0 = y[2*n_surf]
thetf0 = thet[2*n_surf]
mf0 = tan(thetf0)

# final secondary ray parameters
zf1 = z[2*n_surf+1]
yf1 = y[2*n_surf+1]
thetf1 = thet[2*n_surf+1]
mf1 = tan(thetf1)

# compute intercepts
screen_intercept0 = yf0 + mf0*(screen_pos - zf0)
z_intercept0 = zf0 - yf0/mf0
screen_intercept1 = yf1 + mf1*(screen_pos - zf1)
z_intercept1 = zf1 - yf1/mf1

print("+y ray intercepts screen at (%.3f, %f); z-axis at (%f, 0.0)" % \
                (screen_pos,screen_intercept0, z_intercept0))
print("-y ray intercepts screen at (%.3f, %f); z-axis at (%f, 0.0)" % \
                (screen_pos,screen_intercept1, z_intercept1))

# calculate ray intersection position
zint = (yf1 - yf0 + mf0*zf0 - mf1*zf1)/(mf0-mf1)
yint = yf0 + mf0*(zint - zf0)

print("Rays intersect at (z,y) = (%f, %f)" % (zint,yint))

