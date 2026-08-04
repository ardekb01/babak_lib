#include <cstddef>
#include <cfloat>
#include <cmath>
#include "babak_lib.h"
#include "interpolator.h"

// Compute the normalized cross-correlation between the image
// and its reflection across the plane
//
//      Ax + By + Cz = 1
//
// using trilinear interpolation.
float symm_objective_func(const short *image, const DIM &dim, float A, float B, float C)
{
   if( image == nullptr )
      return 0.0f;

   if(dim.nx <= 0 ||
      dim.ny <= 0 ||
      dim.nz <= 0 ||
      dim.dx <= 0.0f ||
      dim.dy <= 0.0f ||
      dim.dz <= 0.0f)
   {
      return 0.0f;
   }

   const double planeNorm2 =
      static_cast<double>(A) * A +
      static_cast<double>(B) * B +
      static_cast<double>(C) * C;

   if(planeNorm2 < DBL_EPSILON)
   {
      return 0.0f;
   }

   const double invNorm = 1.0 / std::sqrt(planeNorm2);

   const float a = static_cast<float>(A * invNorm);
   const float b = static_cast<float>(B * invNorm);
   const float c = static_cast<float>(C * invNorm);
   const float d = static_cast<float>(invNorm);

   int i, j, k;
   float x, y, z;
   float dp;
   float planeK; 
   float planeJK; 
   float planeValue;
   size_t q = 0;
   size_t N = 0;

   float g;

   double Sfg = 0.0;
   double Sg = 0.0;
   double Sgg = 0.0;
   double Sff = 0.0;
   double Sf = 0.0;

   const float ic = (dim.nx - 1.0f) * 0.5f;
   const float jc = (dim.ny - 1.0f) * 0.5f;
   const float kc = (dim.nz - 1.0f) * 0.5f;

   const int np = dim.nx * dim.ny;

   const float a1 = a * dim.dx;
   const float b1 = b * dim.dy;
   const float c1 = c * dim.dz;

   const float invDx = 1.0f / dim.dx;
   const float invDy = 1.0f / dim.dy;
   const float invDz = 1.0f / dim.dz;

   const float reflX = a * invDx;
   const float reflY = b * invDy;
   const float reflZ = c * invDz;

   // inexpensive sanity check
   const size_t nv =
      static_cast<size_t>(dim.nx) *
      static_cast<size_t>(dim.ny) *
      static_cast<size_t>(dim.nz);

   if(nv == 0)
   {
      return 0.0f;
   }

   for(k = 0; k < dim.nz; k++)
   {
      planeK = (k - kc) * c1;

      for(j = 0; j < dim.ny; j++)
      {
         planeJK = planeK + (j - jc) * b1;

         for(i = 0; i < dim.nx; i++)
         {
            const float f = static_cast<float>(image[q]);

            planeValue = planeJK + (i - ic) * a1;

            const float distance = d - planeValue;

            // Evaluate only one side of the plane to avoid
            // counting reflected voxel pairs twice.
            if(f != 0.0f && distance < 0.0f)
            {
               dp = 2.0f * distance;

               x = i + reflX * dp;
               y = j + reflY * dp;
               z = k + reflZ * dp;

               g = linearInterpolator(
                  x,
                  y,
                  z,
                  image,
                  dim.nx,
                  dim.ny,
                  dim.nz,
                  np
               );

               if(g != 0.0f)
               {
                  Sfg += f * g;
                  Sgg += g * g;
                  Sg += g;
                  Sff += f * f;
                  Sf += f;

                  N++;
               }
            }

            q++;
         }
      }
   }

   if(N == 0)
   {
      return 0.0f;
   }

   const double varF = Sff - Sf * Sf / N;
   const double varG = Sgg - Sg * Sg / N;

   if(varF <= DBL_EPSILON || varG <= DBL_EPSILON)
   {
      return 0.0f;
   }

   const double cov = Sfg - Sf * Sg / N;

   float ncc = static_cast<float>(cov / std::sqrt(varF * varG));

   if(ncc > 1.0f)
   {
      ncc = 1.0f;
   }
   else if(ncc < -1.0f)
   {
      ncc = -1.0f;
   }

   return ncc;
}
