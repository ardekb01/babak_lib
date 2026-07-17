#include <stdlib.h>
#include <minmax.h>
#include <babak_lib.h>

void setLowHigh(short *image, int nv, int *low, int *high)
{
   short min = 0;
   short max = 0;

   int *histogram;
   int hsize;        /* histogram size */
   int b;
   int i;

   int nmax;
   int n;

   minmax(image, nv, min, max);

   hsize = max - min + 1;

   histogram = (int *)calloc(hsize, sizeof(int));

   for(i = 0; i < nv; i++)
   {
      b = image[i] - min;

      if(b >= 0 && b < hsize)
         histogram[b]++;
   }

   nmax = (int)(0.0005 * nv);

   n = 0;
   for(i = 0; i < hsize; i++)
   {
      n += histogram[i];

      if(n > nmax)
         break;
   }

   *low = i + min;

   n = 0;
   for(i = 0; i < hsize; i++)
   {
      n += histogram[hsize - 1 - i];

      if(n > nmax)
         break;
   }

   *high = hsize - 1 - i + min;

   free(histogram);
}


void setLowHigh(short *image, int nv, int *low, int *high, float4 percent)
{
   short min = 0;
   short max = 0;

   int *histogram;
   int hsize;        /* histogram size */
   int b;
   int i;

   int nmax;
   int n;

   minmax(image, nv, min, max);

   hsize = max - min + 1;

   histogram = (int *)calloc(hsize, sizeof(int));

   for(i = 0; i < nv; i++)
   {
      b = image[i] - min;

      if(b >= 0 && b < hsize)
         histogram[b]++;
   }

   nmax = (int)(percent * nv / 100.0);

   n = 0;
   for(i = 0; i < hsize; i++)
   {
      n += histogram[i];

      if(n > nmax)
         break;
   }

   *low = i + min;

   n = 0;
   for(i = 0; i < hsize; i++)
   {
      n += histogram[hsize - 1 - i];

      if(n > nmax)
         break;
   }

   *high = hsize - 1 - i + min;

   free(histogram);
}

void setMX(short *image, short *msk, int nv, int &high, float4 alpha)
{
   short min = 0;
   short max = 0;

   int *histogram;
   float4 nmax;

   /*
      Set everything outside the mask equal zero
      and remove negative values if any from the image.
   */
   for(int i = 0; i < nv; i++)
   {
      if(msk[i] == 0)
         image[i] = 0;

      if(image[i] < 0)
         image[i] = 0;
   }


   /*
      Find the maximum of image[] within the mask.
      Minimum is always zero.
   */
   minmax(image, nv, min, max);


   /*
      Allocate memory for histogram.
      Histogram size allows indices:
      histogram[0], histogram[1], ..., histogram[max].
   */
   histogram = (int *)calloc(max + 1, sizeof(int));


   /*
      Fill the histogram.
   */
   {
      int b;

      for(int i = 0; i < nv; i++)
      {
         b = image[i];

         /* Extra precaution to ensure index is not out of range. */
         if(b >= 0 && b <= max)
            histogram[b]++;
      }
   }


   /*
      (nv - histogram[0]) should equal mask size.
   */
   int msksize = nv - histogram[0];

   nmax = alpha * msksize;


   /*
      Find high intensity threshold.
   */
   {
      int i;
      int n;

      n = 0;

      for(i = 0; i <= max; i++)
      {
         n += histogram[max - i];

         if(n > nmax)
            break;
      }

      high = max - i;
   }

   free(histogram);
}
