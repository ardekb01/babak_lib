#include <babak_lib.h>

void compute_cm(short *image, int nx, int ny, int nz, float dx, float dy, float dz, float *x, float *y, float *z)
{
   double dum;
   int q;

   *x=0.0; *y=0.0; *z=0.0;

   dum=0.0;
   q=0;
   for(int k=0;k<nz;k++)
   for(int j=0;j<ny;j++)
   for(int i=0;i<nx;i++)
   {
      dum += image[q];
      *x += i*image[q];
      *y += j*image[q];
      *z += k*image[q];
      q++;
   }

   *x = *x/dum;
   *y = *y/dum;
   *z = *z/dum;

   // transfer origin to the volume center 
   *x -= (nx-1.0)/2.0;
   *y -= (ny-1.0)/2.0;
   *z -= (nz-1.0)/2.0;

   // convert to mm 
   *x *= dx;
   *y *= dy;
   *z *= dz;
}

void compute_cm(short *image, DIM dim, float *cm)
{
   double dum;
   int q;

   cm[0]=0.0; cm[1]=0.0; cm[2]=0.0;

   dum=0.0;
   q=0;
   for(int k=0;k<dim.nz;k++)
   for(int j=0;j<dim.ny;j++)
   for(int i=0;i<dim.nx;i++)
   {
      dum += image[q];
      cm[0] += i*image[q];
      cm[1] += j*image[q];
      cm[2] += k*image[q];
      q++;
   }

   cm[0] = cm[0]/dum;
   cm[1] = cm[1]/dum;
   cm[2] = cm[2]/dum;

   // transfer origin to the volume center 
   cm[0] -= (dim.nx-1.0)/2.0;
   cm[1] -= (dim.ny-1.0)/2.0;
   cm[2] -= (dim.nz-1.0)/2.0;

   // convert to mm 
   cm[0] *= dim.dx;
   cm[1] *= dim.dy;
   cm[2] *= dim.dz;
}

// I represnts a 3x3 symmetric matrix represetned as a 6x1 vector the usual way,
// that is, I[0]=[0][0] I[1]=[1][1] I[2]=[2][2] I[3]=[1][0] I[4]=[2][1] I[5]=[2][0]
void compute_cov(short *image, DIM dim, float *cm, double *I)
{
  float x,y,z; 
  double total_mass;
  int q;

  compute_cm(image, dim, cm);

  for(int i=0; i<6; i++) I[i]=0.0;

  total_mass=0.0;
  q=0;
  for(int k=0;k<dim.nz;k++)
  for(int j=0;j<dim.ny;j++)
  for(int i=0;i<dim.nx;i++)
  {
    total_mass += image[q];
    x = (i - (dim.nx-1.0)/2.0 ) * dim.dx - cm[0];
    y = (j - (dim.ny-1.0)/2.0 ) * dim.dy - cm[1];
    z = (k - (dim.nz-1.0)/2.0 ) * dim.dz - cm[2];

    I[0] += x*x*image[q];
    I[1] += y*y*image[q];
    I[2] += z*z*image[q];
    I[3] += x*y*image[q];
    I[4] += y*z*image[q];
    I[5] += x*z*image[q];

    q++;
  }

  for(int i=0; i<6; i++) I[i] /= total_mass;
}
