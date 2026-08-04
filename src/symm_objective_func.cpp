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
float symm_objective_func(short *image, const DIM &dim, float A, float B, float C)
{
   float a,b,c,d;
   float dum = 0.0f;

   dum = A * A + B * B + C * C;

   if( dum < FLT_EPSILON)
   {
      return(0.0f);
   }

   dum = sqrtf(dum);
   a = A/dum; b = B/dum; c = C/dum;
   d = 1/dum;

   int i, j, k;
   float x, y, z;
   float dp;
   float dum1, dum2, dum3;
   int q;
   int N;

   float f, g;

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

   N = q = 0;

   const float invDx = 1.0f / dim.dx;
   const float invDy = 1.0f / dim.dy;
   const float invDz = 1.0f / dim.dz;

   for(k = 0; k < dim.nz; k++)
   {
      dum1 = (k - kc) * c1;

      for(j = 0; j < dim.ny; j++)
      {
         dum2 = dum1 + (j - jc) * b1;

         for(i = 0; i < dim.nx; i++)
         {
            f = image[q];

            dum3 = dum2 + (i - ic) * a1;

            /* if(f > thresh && (d - dum3) < 0) */
            if(f != 0.0f && (d - dum3) < 0.0f)
            {
               dp = 2.0f * (d - dum3);

               x = i + a * dp * invDx;
               y = j + b * dp * invDy;
               z = k + c * dp * invDz;

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

               /* if(g > thresh) */
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

   return static_cast<float>(cov / std::sqrt(varF * varG));
}

// Compute the normalized cross-correlation between the image
// and its reflection across the plane
//
//      ax + by + cz = d
//
// using trilinear interpolation.
float symm_objective_func(short *image, const DIM &dim, float a, float b, float c, float d)
{
   int i, j, k;
   float x, y, z;
   float dp;
   float dum1, dum2, dum3;
   int q;
   int N;

   float f, g;

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

   N = q = 0;

   const float invDx = 1.0f / dim.dx;
   const float invDy = 1.0f / dim.dy;
   const float invDz = 1.0f / dim.dz;

   for(k = 0; k < dim.nz; k++)
   {
      dum1 = (k - kc) * c1;

      for(j = 0; j < dim.ny; j++)
      {
         dum2 = dum1 + (j - jc) * b1;

         for(i = 0; i < dim.nx; i++)
         {
            f = image[q];

            dum3 = dum2 + (i - ic) * a1;

            /* if(f > thresh && (d - dum3) < 0) */
            if(f != 0.0f && (d - dum3) < 0.0f)
            {
               dp = 2.0f * (d - dum3);

               x = i + a * dp * invDx;
               y = j + b * dp * invDy;
               z = k + c * dp * invDz;

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

               /* if(g > thresh) */
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

   return static_cast<float>(cov / std::sqrt(varF * varG));
}
