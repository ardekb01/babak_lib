#include <stdlib.h>
#include <minmax.h>
#include <babak_lib.h>

void setLowHigh(int2 *image, int4 nv, int4 *low, int4 *high)
{
   int2 min=0, max=0;
   int4 *histogram;
   int4 hsize;			/* histogram size */
   int4 b;
   int4 i;

   int4 nmax;
   int4 n;

   minmax(image,nv,min,max);

   hsize = max-min+1;
   
   histogram=(int4 *)calloc(hsize,sizeof(int4));

   for(i=0;i<nv;i++)
   {
      b = image[i] - min;

      if(b>=0 && b<hsize) histogram[ b ]++;
   }

   nmax = (int4)(0.0005 * nv);

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[i];
      if(n>nmax) break;
   }
   *low=i+min; 

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[hsize-1-i];
      if(n>nmax) break;
   }
   *high=hsize-1-i+min; 

   free(histogram);
}

void setLowHigh(int2 *image, int4 nv, int4 *low, int4 *high, float4 percent)
{
   int2 min=0, max=0;
   int4 *histogram;
   int4 hsize;			/* histogram size */
   int4 b;
   int4 i;

   int4 nmax;
   int4 n;

   minmax(image,nv,min,max);

   hsize = max-min+1;

   histogram=(int4 *)calloc(hsize,sizeof(int4));

   for(int4 i=0; i<nv; i++)
   {
      b = image[i] - min;

      if(b>=0 && b<hsize) histogram[ b ]++;
   }

   nmax = (int4)( percent * nv/100.0);

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[i];
      if(n>nmax) break;
   }
   *low=i+min; 

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[hsize-1-i];
      if(n>nmax) break;
   }
   *high=hsize-1-i+min; 

   free(histogram);
}

void setMX(int2 *image, int2 *msk, int4 nv, int4 &high, float4 alpha)
{
   int2 min=0, max=0;
   int4 *histogram;
   float4 nmax;

   // Set everything outside the msk equal zero
   // and remove negative values if any from the image
   for(int4 i=0; i<nv; i++) 
   {
      if(msk[i]==0) image[i]=0;
      if(image[i]<0) image[i]=0;
   }

   // find the maximum of image[] within the mask
   // minimum is always zero
   minmax(image,nv,min,max);

   // allocate memory for histogram
   // set histogram size so we can have histogram[0], histogram[1], ..., histogram[max]
   histogram=(int4 *)calloc(max+1,sizeof(int4));

   // fill the histogram
   {
      int4 b;
      for(int4 i=0; i<nv; i++)
      {
         b = image[i];

         // extra precaution to ensure index is not out of range 
         if(b>=0 && b<=max ) histogram[ b ]++;
      }
   }

   // (nv-histogram[0]) should equal msk size
   int4 msksize=nv-histogram[0];
   
   nmax = alpha*msksize;

   // find *high
   {
      int4 i;
      int4 n;
      n=0;
      for(i=0; i<=max; i++)
      {
         n += histogram[max-i];
         if(n>nmax) break;
      }

      high=max-i; 
   }

   free(histogram);
}
