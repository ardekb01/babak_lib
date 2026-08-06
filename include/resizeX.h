#ifndef _resize_h
#include <cstddef>
#include "bbk.h"

template<class TYPE> float *resizeX(TYPE *image1, int nx1, int ny1, int nz1, float dx1, int nx2, float dx2)
{
   int n;
   int np1,np2,nv2;
   int i0;
   float d,dl;
   float sd=0.0;
   float *h;
   float *image2;
   float xc1,xc2;
   TYPE *x;

   xc1 = dx1*(nx1-1.0)/2.0;
   xc2 = dx2*(nx2-1.0)/2.0;

   np1=nx1*ny1;
   np2=nx2*ny1;
   nv2=np2*nz1;

   image2=(float *)calloc(nv2,sizeof(float));
   if(image2==NULL) return(NULL);

   if(dx1 < dx2)
   {
      sd=(float)sqrt( (0.5/log(2.0)) * ( dx2*dx2 - dx1*dx1 ) );

      sd /= dx1;
   }
   else sd=0.0;

   h = gaussian_kernel(sd,&n);

   for(int k=0;k<nz1;k++)
   for(int j=0;j<ny1;j++)
   {
      x=image1 + np1*k + nx1*j;
      for(int i=0;i<nx2;i++)
      {
         d = (i*dx2 - xc2 + xc1)/dx1;
         i0=(int)d;
         dl=d - i0;
   
         if(dl==0.0)
            image2[np2*k+nx2*j+i]=conv_pnt_sk(x,nx1,h,n,i0);
         else
            image2[np2*k+nx2*j+i]=(1.0-dl)*conv_pnt_sk(x,nx1,h,n,i0) + dl*conv_pnt_sk(x,nx1,h,n,i0+1);
      }
   }

   free(h);

   return(image2);
}

template<class TYPE> float *resizeX(TYPE *image1, int nx1, int ny1, float dx1, int nx2, float dx2)
{
   int n;
   int np2;
   int i0;
   float d,dl;
   float sd=0.0;
   float *h;
   float *image2;
   float xc1,xc2;
   TYPE *x;

   xc1 = dx1*(nx1-1.0)/2.0;
   xc2 = dx2*(nx2-1.0)/2.0;

   np2=nx2*ny1;

   image2=(float *)calloc(np2,sizeof(float));
   if(image2==NULL) return(NULL);

   if(dx1 < dx2)
   {
      sd=(float)sqrt( (0.5/log(2.0)) * ( dx2*dx2 - dx1*dx1 ) );

      sd /= dx1;
   }
   else sd=0.0;
      
   h = gaussian_kernel(sd,&n);

   for(int j=0;j<ny1;j++)
   {
      x = image1 + nx1*j;

      for(int i=0;i<nx2;i++)
      {
         d = (i*dx2 - xc2 + xc1)/dx1;
         i0=(int)d;
         dl=d - i0;

         if(dl==0.0)
            image2[nx2*j+i]=conv_pnt_sk(x,nx1,h,n,i0);
         else
            image2[nx2*j+i]=(1.0-dl)*conv_pnt_sk(x,nx1,h,n,i0) + dl*conv_pnt_sk(x,nx1,h,n,i0+1);
      }
   }

   free(h);

   return(image2);
}

#define _resize_h
#endif
