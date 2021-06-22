#include <stdio.h>			// for printf, sscanf
#include <stdlib.h>			// for exit()
#include <math.h>			// for trig, sqrt

void ray_step(double *in_ray, double *surface, double *out_ray);

int main(int argc, char* argv[])
{

  double z[40],y[40],thet[40],n[20],z_vert[20],R[20],K[20];
  double in_ray[3],out_ray[3],surf_params[5];
  double startz,starty,slope,offset,current_z;
  double screen_pos=0.0,screen_intercept0,screen_intercept1;
  double z_intercept0,z_intercept1;
  double zf0,yf0,thetf0,mf0,zf1,yf1,thetf1,mf1,zint,yint;
  int n_surf,i;
  char filename[40];
  FILE *lens_file;

  if (argc > 5)			// must have at least these five
  {				// command line arguments
    sscanf(argv[1],"%s",filename);
    sscanf(argv[2],"%lf",&startz);
    sscanf(argv[3],"%lf",&starty);
    sscanf(argv[4],"%lf",&slope);
    sscanf(argv[5],"%lf",&offset);
  }
  else
  {
    printf("Must supply lens_file_name z0, y0, slope, offset arguments\n");
    exit(1);
  }
  if (argc > 6)		// optionally, put a screen somewhere
  {
    sscanf(argv[6],"%lf",&screen_pos);
  }

  // define +y ray to be one offset above start position
  z[0] = startz;
  y[0] = starty + offset;
  thet[0] = atan(slope);	// input slope --> angle in radians

  // define -y ray to be one offset below start position; same angle
  z[1] = startz;
  y[1] =  starty - offset;
  thet[1] = thet[0];

  lens_file = fopen(filename,"r");	// grab lens surface parameters
  fscanf(lens_file,"%d",&n_surf);	// number of surfaces (first line)
  fscanf(lens_file,"%lf",&n[0]);	// initial refr. index (2nd line)
  current_z = 0.0;
  for (i = 0; i < n_surf; i++)		// and n_surf additional lines...
  {
    // read in lens file and verify results
    fscanf(lens_file,"%lf %lf %lf %lf",&n[i+1],&z_vert[i+1],&R[i+1],&K[i+1]);
    current_z += z_vert[i+1];
    printf("Surface %d has n = %f, z_vert = %f, radius = %g, K = %f\n",
           i+1,n[i+1],current_z,R[i+1],K[i+1]);
  }
  close(lens_file);

  printf("Ray 1+ has z = %f; y = %f; thet = %f\n",z[0],y[0],thet[0]);
  printf("Ray 1- has z = %f; y = %f; thet = %f\n",z[1],y[1],thet[1]);

  current_z = 0.0;
  for (i = 0; i < n_surf; i++)	// now propagate surface-at-a-time
  {				// begin ray propagation for loop
    // populate surface parameters array
    current_z += z_vert[i+1];
    surf_params[0] = n[i];
    surf_params[1] = current_z;
    surf_params[2] = R[i+1];
    surf_params[3] = K[i+1];
    surf_params[4] = n[i+1];

    // populate in_ray array for +y ray
    in_ray[0] = z[2*i];
    in_ray[1] = y[2*i];
    in_ray[2] = thet[2*i];

    // carry out +y calculation
    ray_step(in_ray,surf_params,out_ray);

    // stow out_ray into approp. arrays
    z[2*i+2] = out_ray[0];
    y[2*i+2] = out_ray[1];
    thet[2*i+2] = out_ray[2];

    // populate in_ray array for -y ray
    in_ray[0] = z[2*i+1];
    in_ray[1] = y[2*i+1];
    in_ray[2] = thet[2*i+1];

    // carry out -y calculation
    ray_step(in_ray,surf_params,out_ray);

    // stow out_ray into approp. arrays
    z[2*i+3] = out_ray[0];
    y[2*i+3] = out_ray[1];
    thet[2*i+3] = out_ray[2];

    printf("Ray %d+ has z = %f; y = %f; thet = %f\n",
                i+2,z[2*i+2],y[2*i+2],thet[2*i+2]);
    printf("Ray %d- has z = %f; y = %f; thet = %f\n",
                i+2,z[2*i+3],y[2*i+3],thet[2*i+3]);
  }	// end ray propagation for loop

  // final +y ray parameters
  zf0 = z[2*n_surf];
  yf0 = y[2*n_surf];
  thetf0 = thet[2*n_surf];
  mf0 = tan(thetf0);

  // final secondary ray parameters
  zf1 = z[2*n_surf+1];
  yf1 = y[2*n_surf+1];
  thetf1 = thet[2*n_surf+1];
  mf1 = tan(thetf1);

  // compute intercepts
  screen_intercept0 = yf0 + mf0*(screen_pos - zf0);
  z_intercept0 = zf0 - yf0/mf0;
  screen_intercept1 = yf1 + mf1*(screen_pos - zf1);
  z_intercept1 = zf1 - yf1/mf1;

  printf("+y ray intercepts screen at (%.3f, %f); z-axis at (%f, 0.0)\n",
          screen_pos,screen_intercept0, z_intercept0);
  printf("-y ray intercepts screen at (%.3f, %f); z-axis at (%f, 0.0)\n",
          screen_pos,screen_intercept1, z_intercept1);

  // calculate ray intersection position
  zint = (yf1 - yf0 + mf0*zf0 - mf1*zf1)/(mf0-mf1);
  yint = yf0 + mf0*(zint - zf0);

  printf("Rays intersect at (z,y) = (%f, %f)\n",zint,yint);

}

void ray_step(double *in_ray, double *surface, double *out_ray)
{
  double z0,y0,thet0,n0,z_v,R,K,Rp,n_new,Rsign,Ksign;
  double m,denom,b,c,thet_norm,thet_in,thet_out;
  double z_new,y_new,thet_new;

  // rename in_ray components for clarity
  z0 = in_ray[0];
  y0 = in_ray[1];
  thet0 = in_ray[2];

  // rename surface params for clarity
  n0 = surface[0];
  z_v = surface[1];
  R = surface[2];
  K = surface[3];
  n_new = surface[4];

  // establish sign variable for sign flips
  Rsign = R/fabs(R);
  Ksign = 1.0;
  if (K < -1.0) Ksign *= -1.0;

  // compute solution to quadratic formula
  m = tan(thet0);
  if (K == -1.0 && m == 0.0)	// would make denom = 0
  {
    z_new = z_v + y0*y0/(2.0*R);	// parabolic case
  } else {
    denom = 1.0 + K + m*m;
    b = (m*y0 - m*m*z0 - R - (K+1)*z_v)/denom;
    c = (y0*y0 + m*m*z0*z0 - 2.0*m*z0*y0 + 2.0*R*z_v + (K+1)*z_v*z_v)/denom;
    z_new = -b - Rsign*Ksign*sqrt(b*b - c);
  }
  if (fabs(R) > 1.0e10) z_new = z_v;	// planar case

  y_new = y0 + m*(z_new - z0);

  thet_norm = atan(-y_new*Rsign/sqrt(R*R - (K+1)*y_new*y_new));
  thet_in = thet0 - thet_norm;
  thet_out = asin(n0*sin(thet_in)/n_new);
  thet_new = thet_out + thet_norm;
  //printf("%f %f %f %f\n",z_v,z_new,y_new,thet_norm);

  if (z_new < z0) printf("WARNING: ray jumped backwards!\n");

  out_ray[0] = z_new;
  out_ray[1] = y_new;
  out_ray[2] = thet_new;
  
}
