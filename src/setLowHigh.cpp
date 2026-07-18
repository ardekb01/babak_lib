#include "minmax.h"

#include <stdlib.h>

// Computes lower and upper intensity thresholds by discarding
// approximately 'percent' of the voxels from each tail of the intensity histogram.
//
// For example: after calling this function, low and high may be used as 
// follows:
//
//    if(image[i] <= low || image[i] >= high)
//       image[i] = 0;
// This function does not modify the image.
void setLowHigh(const short *image, int nv, int &low, int &high, float percent)
{
   short min = 0;
   short max = 0;

   int hsize;        /* histogram size */
   int b;
   int i;

   int nmax;
   int n;

   // Clamp percent to a valid range.
   // Restricting 'percent' to [0, 50] is sensible and prevents misuse.
   if(percent < 0.0f)
      percent = 0.0f;
   else if(percent > 50.0f)
      percent = 50.0f;

   // Validate input arguments.
   if(image == nullptr || nv <= 0)
   {
      low = 0;
      high = 0;
      return;
   }

   // Find the minimum and maximum voxel intensities.
   // A signed 16-bit short ranges from -32768 to 32767.
   minmax(image, nv, min, max);

   // Maximum histogram size for a signed 16-bit image is 32767 - (-32768) + 1 = 65536 bins
   hsize = max - min + 1;

   int *histogram = (int *)calloc(hsize, sizeof(int));

   // guard against memory allocation failure
   if(histogram == nullptr)
   {
      low = min;
      high = max;
      return;
   }

   for(i = 0; i < nv; i++)
   {
      b = image[i] - min;

      // This should always be true, but keep the check as a safeguard.
      if(b >= 0 && b < hsize)
         histogram[b]++;
   }

   nmax = (int)(percent * nv / 100.0f + 0.5f);

   n = 0;
   for(i = 0; i < hsize; i++)
   {
      n += histogram[i];

      if(n > nmax)
         break;
   }

   low = i + min;

   n = 0;
   for(i = 0; i < hsize; i++)
   {
      n += histogram[hsize - 1 - i];

      if(n > nmax)
         break;
   }

   high = hsize - 1 - i + min;

   free(histogram);
}
