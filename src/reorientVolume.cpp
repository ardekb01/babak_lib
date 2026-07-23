#include <stdlib.h>

#include "babak_lib.h"
#include "nifti1.h"
#include "PILtransform.h"
#include "getNiftiImageOrientation.h"

// The upper-left 3x3 submatrix of T is a signed permutation matrix.
// Callers must release the returned image with free()
short *reorientVolume(short *input_image,
                       int nx1,
                       int ny1,
                       int nz1,
                       float dx1,
                       float dy1,
                       float dz1,
                       float *T,
                       int &nx2,
                       int &ny2,
                       int &nz2,
                       float &dx2,
                       float &dy2,
                       float &dz2)
{
   // The output dimensions are initialized.
   nx2 = 0;
   ny2 = 0;
   nz2 = 0;

   dx2 = 0.0f;
   dy2 = 0.0f;
   dz2 = 0.0f;

   // Validate input variables.
   if (input_image == nullptr || T == nullptr)
   {
      return nullptr;
   }

   // Dimension and voxel-size validation.
   if (nx1 <= 0 || ny1 <= 0 || nz1 <= 0 ||
       dx1 <= 0.0f || dy1 <= 0.0f || dz1 <= 0.0f)
   {
      return nullptr;
   }

   // Validate that T is a signed permutation matrix.
   for (int i = 0; i < 3; i++)
   {
      int count = 0;

      for (int j = 0; j < 3; j++)
      {
         float value = T[4 * i + j];

         if (value != 0.0f)
         {
            if (value != 1.0f && value != -1.0f)
            {
               return nullptr;
            }

            count++;
         }
      }

      if (count != 1)
      {
         return nullptr;
      }
   }

   for (int j = 0; j < 3; j++)
   {
      int count = 0;

      for (int i = 0; i < 3; i++)
      {
         if (T[4 * i + j] != 0.0f)
         {
            count++;
         }
      }

      if (count != 1)
      {
         return nullptr;
      }
   }

   // Set nx2, ny2, nz2, dx2, dy2, and dz2.
   if (T[0] != 0.0f)
   {
      nx2 = nx1;
      dx2 = dx1;
   }
   else if (T[4] != 0.0f)
   {
      ny2 = nx1;
      dy2 = dx1;
   }
   else if (T[8] != 0.0f)
   {
      nz2 = nx1;
      dz2 = dx1;
   }

   if (T[1] != 0.0f)
   {
      nx2 = ny1;
      dx2 = dy1;
   }
   else if (T[5] != 0.0f)
   {
      ny2 = ny1;
      dy2 = dy1;
   }
   else if (T[9] != 0.0f)
   {
      nz2 = ny1;
      dz2 = dy1;
   }

   if (T[2] != 0.0f)
   {
      nx2 = nz1;
      dx2 = dz1;
   }
   else if (T[6] != 0.0f)
   {
      ny2 = nz1;
      dy2 = dz1;
   }
   else if (T[10] != 0.0f)
   {
      nz2 = nz1;
      dz2 = dz1;
   }

   // Allocate memory for the output image.
   size_t nvoxels = (size_t)nx2 * ny2 * nz2;

   short *output_image = (short *)malloc(nvoxels * sizeof(short));

   if (output_image == nullptr)
   {
      return nullptr;
   }

   int xc;
   int yc;
   int zc;

   // If necessary, reverse an axis.
   if (T[0] < 0.0f || T[4] < 0.0f || T[8] < 0.0f)
   {
      xc = nx1 - 1;
   }
   else
   {
      xc = 0;
   }

   if (T[1] < 0.0f || T[5] < 0.0f || T[9] < 0.0f)
   {
      yc = ny1 - 1;
   }
   else
   {
      yc = 0;
   }

   if (T[2] < 0.0f || T[6] < 0.0f || T[10] < 0.0f)
   {
      zc = nz1 - 1;
   }
   else
   {
      zc = 0;
   }

   size_t q = 0;
   size_t np1 = (size_t)nx1 * ny1;

   int i1_a;
   int j1_a;
   int k1_a;

   int i1_b;
   int j1_b;
   int k1_b;

   int i1;
   int j1;
   int k1;

   size_t input_index;

   for (int k2 = 0; k2 < nz2; k2++)
   {
      i1_a = (int)T[8] * k2 + xc;
      j1_a = (int)T[9] * k2 + yc;
      k1_a = (int)T[10] * k2 + zc;

      for (int j2 = 0; j2 < ny2; j2++)
      {
         i1_b = (int)T[4] * j2 + i1_a;
         j1_b = (int)T[5] * j2 + j1_a;
         k1_b = (int)T[6] * j2 + k1_a;

         for (int i2 = 0; i2 < nx2; i2++)
         {
            i1 = (int)T[0] * i2 + i1_b;
            j1 = (int)T[1] * i2 + j1_b;
            k1 = (int)T[2] * i2 + k1_b;

            input_index = (size_t)k1 * np1 +
                          (size_t)j1 * nx1 +
                          (size_t)i1;

            output_image[q++] = input_image[input_index];
         }
      }
   }

   return output_image;
}
