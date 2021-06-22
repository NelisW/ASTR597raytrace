from math import *
import numpy
from scipy.optimize import leastsq

def vmag(vect):
    return sqrt(numpy.dot(vect,vect))

def vnorm(vect):
    return vect/vmag(vect)

def ray_step_3d(ray_pos, ray_direc, surface):

    # rename ray_pos components for clarity
    x0 = ray_pos[0]
    y0 = ray_pos[1]
    z0 = ray_pos[2]
    init_pos = numpy.array([x0,y0,z0],dtype='d')

    # rename ray_direc components for clarity, and normalize
    kx0 = ray_direc[0]
    ky0 = ray_direc[1]
    kz0 = ray_direc[2]
    k0 = vnorm(numpy.array([kx0,ky0,kz0],dtype='d'))

    # rename surface params for clarity
    n0 = surface[0]
    z_v = surface[1]
    R = surface[2]
    K = surface[3]
    n_new = surface[4]

    # establish sign flippy dealies
    Rsign = R/fabs(R)
    Ksign = 1.0
    if (K < -1.0):
        Ksign = -1.0
    direc_sign = kz0/fabs(kz0)

    # solve intersection
    denom = kx0*kx0 + ky0*ky0
    if (K == -1.0 and denom == 0.0):		# parabolic case, straight in
        s = (x0*x0 + y0*y0 - 2*R*(z0 - z_v))/(2*R*kz0)
    else:
        denom = kx0*kx0 + ky0*ky0 + (K+1)*kz0*kz0
        b = (x0*kx0 + y0*ky0 + ((K+1)*(z0-z_v) - R)*kz0)/denom
        c = (x0*x0 + y0*y0 + (K+1)*(z0*z0-2*z0*z_v+z_v*z_v) - 2*R*(z0-z_v))/denom
        s = -b - direc_sign*Rsign*Ksign*sqrt(b*b - c)

    if (fabs(R) > 1.0e10):	# assume exactly planar
        s = (z_v - z0)/kz0

    if (s < 0.0):
        #print "WARNING: ray jumped backwards!"
        backwards = True
    else:
        backwards = False
    
    # ray/surface intersection position
    new_pos = init_pos + s*k0
    xi = new_pos[0]
    yi = new_pos[1]
    zi = new_pos[2]
    ri2 = xi*xi + yi*yi
    
    # surface normal
    nhat = vnorm(numpy.array([xi,yi,-Rsign*sqrt(R*R - (K+1)*ri2)],dtype='d'))
    
    if (fabs(R) > 1.0e10):	# planar case
        nhat = numpy.array([0.0,0.0,1.0],dtype='d')

    # establish delta-theta
    ndotk = numpy.dot(nhat,k0)
    sin_thet_in = sqrt(1.0 - ndotk*ndotk)
    thet_in = asin(sin_thet_in)
    thet_out = asin(n0*sin_thet_in/n_new)
    dthet = thet_out - thet_in

    # establish l-vector perp to k-vect and in plane of incidence
    lvect = nhat - ndotk*k0
    if (vmag(lvect) > 1.0e-6):
        lhat = vnorm(lvect)
    else:
        lhat = numpy.array([0.0,0.0,0.0],dtype='d')

    # get sign flips
    dnsign = n0*n_new/fabs(n0*n_new)
    knsign = -ndotk/fabs(ndotk)

    k_new = dnsign*(cos(dthet)*k0 + knsign*sin(dthet)*lhat)

    return (new_pos,k_new,backwards)

def    ray_step_2d(in_ray, surface):

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
    #print "%f %f %f %f" % (z_v,z_new,y_new,thet_norm)

    if (z_new < z0):
        print "WARNING: ray jumped backwards!"

    out_ray[0] = z_new
    out_ray[1] = y_new
    out_ray[2] = thet_new

    return out_ray

def ray_bundle(ctr_pos,Rmax,approx_n):

    """
    Set up polar grid of aproximately equal-spaced rays
    Someday tilt plane according to central direction?
    """
    x = []
    y = []
    z = []
    if (approx_n < 5):
        approx_n = 5
    density = approx_n/(pi*pow(Rmax,2))
    spacing = sqrt(1.0/density)
    nring = int(round(Rmax/spacing))
    rad_spacing = Rmax/nring
    for ring in range(nring+1):
        rad = ring*rad_spacing
        nray = int(round(2*pi*rad/rad_spacing))
        if (nray == 0):
            nray = 1
        theta = numpy.arange(0.0,2*pi,2*pi/nray)
        for ray in range(nray):
            x.append(ctr_pos[0] + rad*cos(theta[ray]))
            y.append(ctr_pos[1] + rad*sin(theta[ray]))
            z.append(ctr_pos[2])
    Nray = len(x)
    return Nray,x,y,z

def rms_blur(p,state):
    """
         p[0] is z position of screen
         state is [[x], [y], [z], [kx], [ky], [kz]]
    """

    Nray = len(state[0])

    x = state[0]
    y = state[1]
    z = state[2]
    kx = state[3]
    ky = state[4]
    kz = state[5]

    sumx = 0.0
    sumy = 0.0
    sumx2 = 0.0
    sumy2 = 0.0
    for i in range(Nray):
        s = (p[0] - z[i])/kz[i]
        xs = x[i] + s*kx[i]
        ys = y[i] + s*ky[i]
        sumx += xs
        sumy += ys
        sumx2 += xs*xs
        sumy2 += ys*ys
    avgx = sumx/Nray
    avgy = sumy/Nray
    rms = sqrt(sumx2/Nray + sumy2/Nray - avgx*avgx - avgy*avgy)

    return rms

def best_focus(x,y,z,kx,ky,kz):
    """
         Given final ray states, find place with best RMS blur
         p[0] --> z position of screen
    """

    state = [x,y,z,kx,ky,kz]
    pest = numpy.array([0.0],dtype='d')

    fit = leastsq(rms_blur,pest.copy(),args=(state),full_output=1)

    return fit

def centroid(screen,state):
    """
         screen is z position of screen
         state is [[x], [y], [z], [kx], [ky], [kz]]
    """

    Nray = len(state[0])

    x = state[0]
    y = state[1]
    z = state[2]
    kx = state[3]
    ky = state[4]
    kz = state[5]

    sumx = 0.0
    sumy = 0.0
    for i in range(Nray):
        s = (screen - z[i])/kz[i]
        xs = x[i] + s*kx[i]
        ys = y[i] + s*ky[i]
        sumx += xs
        sumy += ys
    avgx = sumx/Nray
    avgy = sumy/Nray

    return avgx,avgy

def spot_pat(screen,state):
    """
         screen is z position of screen
         state is [[x], [y], [z], [kx], [ky], [kz]]
    """

    Nray = len(state[0])

    x = state[0]
    y = state[1]
    z = state[2]
    kx = state[3]
    ky = state[4]
    kz = state[5]

    x_out = []
    y_out = []
    for i in range(Nray):
        s = (screen - z[i])/kz[i]
        x_out.append(x[i] + s*kx[i])
        y_out.append(y[i] + s*ky[i])

    return x_out,y_out

def skewness(screen,state):
    """
         screen is z position of screen
         state is [[x], [y], [z], [kx], [ky], [kz]]
    """

    Nray = len(state[0])

    x = state[0]
    y = state[1]
    z = state[2]
    kx = state[3]
    ky = state[4]
    kz = state[5]

    sumx = 0.0
    sumy = 0.0
    sumx2 = 0.0
    sumy2 = 0.0
    sumx3 = 0.0
    sumy3 = 0.0
    for i in range(Nray):
        s = (screen - z[i])/kz[i]
        xs = x[i] + s*kx[i]
        ys = y[i] + s*ky[i]
        sumx += xs
        sumy += ys
        sumx2 += xs*xs
        sumy2 += ys*ys
        sumx3 += xs*xs*xs
        sumy3 += ys*ys*ys
    avgx = sumx/Nray
    avgy = sumy/Nray
    avgx2 = sumx2/Nray
    avgy2 = sumy2/Nray
    skew_x = (sumx3/Nray - 3*avgx2*avgx + 2*pow(avgx,3))/pow(avgx2-avgx*avgx,1.5)
    skew_y = (sumy3/Nray - 3*avgy2*avgy + 2*pow(avgy,3))/pow(avgy2-avgy*avgy,1.5)

    return skew_x,skew_y
