#include <stdlib.h>
#include <stdio.h>
#include "babak_lib.h"
#include "bbk.h"

// NOTE: The standard deviation (SD) in all these functions is in units of pixels
// use smoothY(float *image, int nx, int ny, int nz, float sd) with nz=1 for 2D case

float *smoothY(float *image, int nx, int ny, int nz, float sd)
{
   int n;
   float *h;
   float *image_out;
   float *y;
   int slice_offset;

   image_out=(float *)calloc(nx*ny*nz,sizeof(float));
   if(image_out==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int i=0; i<nx*ny*nz; i++) image_out[i] = image[i]; // bug fix
      return(image_out);
   }

   h = gaussian_kernel(sd, &n);
   if(h==NULL)
   {
      return(NULL);
   }

   y=(float *)calloc(ny,sizeof(float));

   for(int k=0;k<nz;k++)
   {
      slice_offset = nx*ny*k;
      for(int i=0;i<nx;i++)
      {
         for(int j=0;j<ny;j++) 
         {
            y[j]=image[slice_offset + nx*j + i];
         }

         for(int j=0;j<ny;j++)
         {
            image_out[slice_offset + nx*j + i] = conv_pnt_sk(y,ny,h,n,j);
         }
      }
   }

   free(h);
   free(y);

   return(image_out);
}

float *smoothZ(float *image, int nx, int ny, int nz, float sd)
{
   int n;
   int np;
   int row_offset;
   float *h;
   float *z;
   float *image_out;

   np = nx*ny;

   image_out=(float *)calloc(np*nz ,sizeof(float));
   if(image_out==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int i=0; i<nx*ny*nz; i++) image_out[i] = image[i];
      return(image_out);
   }

   h = gaussian_kernel(sd,&n);
   if(h==NULL)
   {
      return(NULL);
   }

   z=(float *)calloc(nz,sizeof(float));

   for(int j=0;j<ny;j++)
   {
      row_offset = nx*j;

      for(int i=0;i<nx;i++)
      {
         for(int k=0;k<nz;k++)
         {
            z[k] = image[np*k + row_offset +i];
         }

         for(int k=0;k<nz;k++)
         {
            image_out[np*k + row_offset +i] = conv_pnt_sk(z,nz,h,n,k);
         }
      }
   }

   free(h);
   free(z);

   return(image_out);
}
