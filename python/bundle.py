#!/usr/bin/env python

from math import *
import sys
import numpy
import pylab
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
surf_params = [0.0]*5

screen_pos=0.0

if (narg > 9):
    filename = sys.argv[1]
    xc = float(sys.argv[2])
    yc = float(sys.argv[3])
    zc = float(sys.argv[4])
    kx0 = float(sys.argv[5])
    ky0 = float(sys.argv[6])
    kz0 = float(sys.argv[7])
    Rmax = float(sys.argv[8])
    nray_target = int(sys.argv[9])
else:
    print "Need lens_file_name x0 y0 z0 kx0 ky0 kz0 Rmax nray arguments"
    sys.exit()

haveOwnScreen = False
if (narg > 10):		# optionally, put a screen somewhere
    haveOwnScreen = True
    screen_pos = float(sys.argv[10])

# initialize ray bundle
Nray,x0,y0,z0 = raytrace.ray_bundle([xc,yc,zc],Rmax,nray_target)
print "Processing %d rays" % Nray

# insert into lists
x.append(x0)
y.append(y0)
z.append(z0)
kx.append([kx0]*Nray)
ky.append([ky0]*Nray)
kz.append([kz0]*Nray)

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
    print "Surface %d has n = %.4f, z_vert = %.3f, radius = %g, K = %.4f" % \
                 (i+1,n[i+1],current_z,R[i+1],K[i+1])

lens_file.close()

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

    # make room for new ray parameters
    x.append([])
    y.append([])
    z.append([])
    kx.append([])
    ky.append([])
    kz.append([])

    nback = 0
    for ray in range(Nray):
        # populate in_ray array for +y ray
        ray_pos = numpy.array([x[i][ray],y[i][ray],z[i][ray]],dtype='d')
        ray_direc = vnorm(numpy.array([kx[i][ray],ky[i][ray],kz[i][ray]],dtype='d'))

        # carry out calculation
        out_pos,out_direc,back = raytrace.ray_step_3d(ray_pos,ray_direc,surf_params)
        if back:
            nback += 1

        # stow out_ray into approp. arrays
        x[i+1].append(out_pos[0])
        y[i+1].append(out_pos[1])
        z[i+1].append(out_pos[2])
        kx[i+1].append(out_direc[0])
        ky[i+1].append(out_direc[1])
        kz[i+1].append(out_direc[2])

    if nback:
        print "%d Rays out of %d jumped backwards at surface %d" % \
                    (nback,Nray,i+1)

if not haveOwnScreen:
    bf_fit = raytrace.best_focus(x[-1],y[-1],z[-1],kx[-1],ky[-1],kz[-1])

    screen_pos = bf_fit[0]

state = [x[-1],y[-1],z[-1],kx[-1],ky[-1],kz[-1]]
rms = raytrace.rms_blur([screen_pos],state)
x_ctr,y_ctr = raytrace.centroid(screen_pos,state)
if haveOwnScreen:
    print "Hits screen at (%f, %f, %f); RMS = %f" % (x_ctr,y_ctr,screen_pos,rms)
else:
    print "Best focus at (%f, %f, %f); RMS = %f" % (x_ctr,y_ctr,screen_pos,rms)
print "Skew in x: %f; skew in y: %f" % raytrace.skewness(screen_pos,state)

x_spot, y_spot = raytrace.spot_pat(screen_pos,state)
ax = pylab.subplot(1,1,1,aspect='equal')
pylab.plot(x_spot,y_spot,'b.',markersize=3)
pylab.savefig('spot.png')
pylab.clf()
