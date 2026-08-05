#include <stdlib.h>
#include <stdio.h>
#include "babak_lib.h"
#include "bbk.h"

// NOTE: The standard deviation (SD) in all these functions is in units of pixels
// use smoothY(float4 *image, int4 nx, int4 ny, int4 nz, float4 sd) with nz=1 for 2D case

float4 *smoothY(float4 *image, int4 nx, int4 ny, int4 nz, float4 sd)
{
   int4 n;
   float4 *h;
   float4 *image_out;
   float4 *y;
   int4 slice_offset;

   image_out=(float4 *)calloc(nx*ny*nz,sizeof(float4));
   if(image_out==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int4 i=0; i<nx*ny*nz; i++) image_out[i] = image[i]; // bug fix
      return(image_out);
   }

   h = gaussian_kernel(sd, &n);
   if(h==NULL)
   {
      return(NULL);
   }

   y=(float4 *)calloc(ny,sizeof(float4));

   for(int4 k=0;k<nz;k++)
   {
      slice_offset = nx*ny*k;
      for(int4 i=0;i<nx;i++)
      {
         for(int4 j=0;j<ny;j++) 
         {
            y[j]=image[slice_offset + nx*j + i];
         }

         for(int4 j=0;j<ny;j++)
         {
            image_out[slice_offset + nx*j + i] = conv_pnt_sk(y,ny,h,n,j);
         }
      }
   }

   free(h);
   free(y);

   return(image_out);
}

float4 *smoothZ(float4 *image, int4 nx, int4 ny, int4 nz, float4 sd)
{
   int4 n;
   int4 np;
   int4 row_offset;
   float4 *h;
   float4 *z;
   float4 *image_out;

   np = nx*ny;

   image_out=(float4 *)calloc(np*nz ,sizeof(float4));
   if(image_out==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int4 i=0; i<nx*ny*nz; i++) image_out[i] = image[i];
      return(image_out);
   }

   h = gaussian_kernel(sd,&n);
   if(h==NULL)
   {
      return(NULL);
   }

   z=(float4 *)calloc(nz,sizeof(float4));

   for(int4 j=0;j<ny;j++)
   {
      row_offset = nx*j;

      for(int4 i=0;i<nx;i++)
      {
         for(int4 k=0;k<nz;k++)
         {
            z[k] = image[np*k + row_offset +i];
         }

         for(int4 k=0;k<nz;k++)
         {
            image_out[np*k + row_offset +i] = conv_pnt_sk(z,nz,h,n,k);
         }
      }
   }

   free(h);
   free(z);

   return(image_out);
}
