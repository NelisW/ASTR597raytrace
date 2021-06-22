#include <stdio.h>			// for printf, sscanf
#include <stdlib.h>			// for exit()
#include <math.h>			// for trig, sqrt

void ray_step(double *in_ray, double *surface, double *out_ray);

int main(int argc, char* argv[])
{

  double z[20],y[20],thet[20],n[20],z_vert[20],R[20],K[20];
  double in_ray[3],out_ray[3],surf_params[5];
  double slope,current_z;
  double screen_pos=0.0,zf,yf,thetf,mf;
  double screen_intercept,z_intercept;
  int n_surf,i;
  char filename[40];
  FILE *lens_file;

  if (argc > 4)			// must have at least these five
  {				// command line arguments
    sscanf(argv[1],"%s",filename);
    sscanf(argv[2],"%lf",&z[0]);
    sscanf(argv[3],"%lf",&y[0]);
    sscanf(argv[4],"%lf",&slope);
  }
  else
  {
    printf("Must supply lens_file_name z0, y0, slope arguments\n");
    exit(1);
  }
  if (argc > 5)		// optionally, put a screen somewhere
  {
    sscanf(argv[5],"%lf",&screen_pos);
  }

  thet[0] = atan(slope);	// input slope --> angle in radians

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

  printf("Ray 1 has z = %f; y = %f; thet = %f\n",z[0],y[0],thet[0]);

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

    // populate in_ray array
    in_ray[0] = z[i];
    in_ray[1] = y[i];
    in_ray[2] = thet[i];

    // carry out calculation
    ray_step(in_ray,surf_params,out_ray);

    // stow out_ray into approp. arrays
    z[i+1] = out_ray[0];
    y[i+1] = out_ray[1];
    thet[i+1] = out_ray[2];

    printf("Ray %d has z = %f; y = %f; thet = %f\n",
                i+2,z[i+1],y[i+1],thet[i+1]);
  }	// end ray propagation for loop

  // final ray parameters
  zf = z[n_surf];
  yf = y[n_surf];
  thetf = thet[n_surf];
  mf = tan(thetf);

  // compute intercepts
  screen_intercept = yf + mf*(screen_pos - zf);
  z_intercept = zf - yf/mf;

  printf("Ray intercepts screen at (%.3f, %f); z-axis at (%f, 0.0)\n",
          screen_pos,screen_intercept, z_intercept);

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
