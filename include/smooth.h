#ifndef _smooth_h
#include <cstddef>
#include "bbk.h"

template<class TYPE> float *smoothX(TYPE *image, int nx, int ny, int nz, float sd)
{
   int n;
   int slice_offset;
   int row_offset;

   float *h;
   float *image_out;
   TYPE *x;

   image_out=(float *)calloc(nx*ny*nz,sizeof(float));
   if(image_out==NULL) return(NULL);

   if(sd<=0.0)
   {
      // for(int i=0; i<nx*ny; i++) image_out[i] = image[i]; // This seemed to be a bug - removed
      for(int i=0; i<nx*ny*nz; i++) image_out[i] = image[i];
      return(image_out);
   }

   h = gaussian_kernel(sd, &n);
   if(h==NULL)
   {
      return(NULL);
   }

   for(int k=0;k<nz;k++)
   {
      slice_offset = nx*ny*k;
      for(int j=0;j<ny;j++)
      {
         row_offset = slice_offset + nx*j;
         for(int i=0;i<nx;i++)
         {
            x = image + row_offset;
            image_out[ row_offset + i ] = conv_pnt_sk(x,nx,h,n,i);
         }
      }
   }

   free(h);

   return(image_out);
}

template<class TYPE> float *smoothXYZ(TYPE *image, int nx, int ny, int nz, float sdx, float sdy, float sdz)
{
   float *imA, *imB, *imC;

   imA=smoothX(image, nx, ny, nz, sdx);
   if(imA==NULL)
   {
      return(NULL);
   }

   imB=smoothY(imA, nx, ny, nz, sdy);
   free(imA);
   if(imB==NULL)
   {
      return(NULL);
   }

   imC=smoothZ(imB, nx, ny, nz, sdz);
   free(imB);

   return(imC);
}

template<class TYPE> float *smoothXY(TYPE *image, int nx, int ny, float sd)
{
   float *imA, *imB;

   imA=smoothX(image, nx, ny, 1, sd);

   imB=smoothY(imA, nx, ny, 1, sd);
   free(imA);

   return(imB);
}
/////////////////////////////////////////////////////////////////////
#define _smooth_h
#endif
