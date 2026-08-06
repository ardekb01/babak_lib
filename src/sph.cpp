#include <babak_lib.h>
#include <sph.h>
#include <math.h>
#include <stdlib.h>

float SPH::dot(float *u)
{
   float dot = 0.0;

   for(int m=0; m<n; m++)
      dot += (v[m]*u[m]);

   return(dot);
}

float SPH::norm()
{
   float norm=0.0;

   for(int m=0; m<n; m++)
      norm += (v[m]*v[m]);

   norm = sqrtf(norm);

   return(norm);
}

float SPH::mean()
{
   float mean=0.0;

   for(int m=0; m<n; m++)
      mean += v[m];

   mean /= n;

   return(mean);
}

void SPH::normalize()
{
   float norm=0.0;

   for(int m=0; m<n; m++)
      norm += (v[m]*v[m]);

   norm = sqrtf(norm);

   for(int m=0; m<n; m++)
      v[m] /= norm;
}

void SPH::zeromean()
{
   float mean=0.0;

   for(int m=0; m<n; m++)
      mean += v[m];

   mean /= n;

   for(int m=0; m<n; m++)
      v[m] -= mean;
}

void SPH::get(float *array)
{
   for(int m=0; m<n; m++)
   {
      array[m] = v[m];
   }
}

void SPH::set(SHORTIM im, int ic, int jc, int kc)
{
   int i0, j0, k0;

   for(int m=0; m<n; m++)
   {
      i0=ic+i[m];
      j0=jc+j[m];
      k0=kc+k[m];

      if(i0>=0 && i0<im.nx && j0>=0 && j0<im.ny && k0>=0 && k0<im.nz)
      {
         v[m] = im.v[ im.np*k0 + im.nx*j0 + i0 ];
      } else
      {
         v[m] = 0.0;
      }
   }
}

void SPH::set(float *im, int nx, int ny, int nz, int ic, int jc, int kc)
{
   int i0, j0, k0;
   int np=nx*ny;

   for(int m=0; m<n; m++)
   {
      i0=ic+i[m];
      j0=jc+j[m];
      k0=kc+k[m];

      if(i0>=0 && i0<nx && j0>=0 && j0<ny && k0>=0 && k0<nz)
      {
         v[m] = im[ np*k0 + nx*j0 + i0 ];
      } else
      {
         v[m] = 0.0;
      }
   }
}

void SPH::reset()
{
   for(int m=0; m<n; m++)
   {
      v[m]=0.0;
   }
}

SPH::~SPH()
{
   delete i;
   delete j;  
   delete k;  
   delete v;  
}

SPH::SPH(int radius)
{
   r=radius;

   int nmax;
   float r2;
   float x2, y2;

   if(r<0) r=0;

   nmax = (2*r+1) * (2*r+1) * (2*r+1);

   i = (int *)calloc(nmax, sizeof(int));
   j = (int *)calloc(nmax, sizeof(int));
   k = (int *)calloc(nmax, sizeof(int));

   r2 = r*r;

   n=0;
   for(int x=-r; x<=r; x++)   
   {
      x2=x*x;
      for(int y=-r; y<=r; y++)   
      {
         y2=y*y;
         for(int z=-r; z<=r; z++)   
         {
            if( (x2+y2+z*z) <= r2 )
            {
               i[n]=x; j[n]=y; k[n]=z;
               n++;
            }
         }
      }
   }

   v = (float *)calloc(n, sizeof(float));
}

