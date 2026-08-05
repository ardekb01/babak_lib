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

   // Validate the input parameters
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

   // norm2 = norm * norm
   const double norm2 =
      static_cast<double>(A) * A +
      static_cast<double>(B) * B +
      static_cast<double>(C) * C;

   if(norm2 <= DBL_EPSILON)
   {
      return 0.0f;
   }

   // The L2 norm of (A, B, C), ||(A,B,C)|| 
   const double norm = std::sqrt(norm2);

   // From the input parameters (A,B,C), we compute the parameters
   // (a,b,c) and d for the equation of the same plane represented
   // as: ax + by + cz = d.  Note ||(a,b,c)|| = 1.0. (a,b,c) is
   // the unit normal to the plane. d is the perpendicular distance
   // from the origin to the plane.
   // Note: (a,b,c).(x,y,z) = d for all points (x,y,z) on the plane.
   // Note: if multiply both sides of the equation Ax + By + Cz = 1
   // by d, we obtained Adx + Bdy + Cdz = d. Therefore, a = A*d,
   // b = B*d and c = C*d, as computed below.
   const float d = static_cast<float>(1.0 / norm);
   const float a = static_cast<float>(A * d);
   const float b = static_cast<float>(B * d);
   const float c = static_cast<float>(C * d);

   // The number of pairs of points on the positive and
   // negative sides of the plane. The normalized cross
   // correlation is computed between these pairs.
   size_t numPairs = 0;

   int i, j, k;
   float dotZ; 
   float dotYZ; 
   float dotXYZ;
   size_t q = 0;

   double Sfg = 0.0;
   double Sg = 0.0;
   double Sgg = 0.0;
   double Sff = 0.0;
   double Sf = 0.0;

   // (ic, jc, kc) are used to take move the origin of the
   // (i, j, k) coordinates system to the center of 
   // the 3D image FOV.  
   // So (i', j', k') = (i-ic, j-jc, k-kc).
   // Multiplication by (dim.dx, dim.dy, dim.dz) 
   // the converts (i', j', k') coordinates to
   // (x, y, z) coordinates, that is:
   // (x, y, z) = (i'*dim.dx, j'*dim.dy, k'*dim.dz)
   // These are the coordinates with respect to which
   // our plane equation Ax + By + Cz = 1 is defined.
   const float ic = (dim.nx - 1.0f) * 0.5f;
   const float jc = (dim.ny - 1.0f) * 0.5f;
   const float kc = (dim.nz - 1.0f) * 0.5f;

   const int np = dim.nx * dim.ny;

   // (a1, b1, c1) and (a2, b2, c2) are slightly modified versions
   // of (a, b, c) defined to save computations inside the loops.
   const float a1 = a * dim.dx;
   const float b1 = b * dim.dy;
   const float c1 = c * dim.dz;

   const float a2 = a / dim.dx;
   const float b2 = b / dim.dy;
   const float c2 = c / dim.dz;

   for(k = 0; k < dim.nz; k++)
   {
      // Computes the dot product: (0,0,z).(a,b,c)
      dotZ = (k - kc) * c1;  

      for(j = 0; j < dim.ny; j++)
      {
         // Computes the dot product: (0,y,z).(a,b,c)
         dotYZ = dotZ + (j - jc) * b1; 

         for(i = 0; i < dim.nx; i++)
         {
            // Computes the dot product: (x,y,z).(a,b,c)
            // This is the signed length of the projection 
            // of point (x,y,z) in the (a,b,c) direction.
            dotXYZ = dotYZ + (i - ic) * a1;

            // Signed distance of the point (x,y,z) to the plane.
            // That is, the signed length of the perpendicular line
            // dropped from point (x,y,z) to the plane. If this
            // distance is positive, then we are on the positive
            // side of the plane.  If it is negative, then we are
            // on the negative side of the plane.  If it is zero, then
            // we are on the plane because then, dotXYZ=d.
            const float distance = d - dotXYZ;

            // f is the value of the image at point (x,y,z).
            const float f = static_cast<float>(image[q]);

            // Evaluate only one side of the plane to avoid
            // counting reflected voxel pairs twice.
            // Only consider voxels with non-zero values. This avoids
            // using background voxels that are not on the head.
            if(f != 0.0f && distance < 0.0f)
            {
               const float offset = 2.0f * distance;

               // (i', j', z') are the coordinates of the reflection
               // (the Householder transformation) of (i, j, k) with respect
               // to our plane.
               // To obtain these, let the reflection of (x, y, z) wrt
               // the plane be (x', y', z').  This is obtained by:
               // (x', y', z') = (x, y, z) + 2(d-dotXYZ)*(a, b, c). Then we
               // convert (x', y', z') back to the (i',j',k') coordinates system,
               // indicated as (ip, jp, kp) below.
               const float ip = i + a2 * offset;
               const float jp = j + b2 * offset;
               const float kp = k + c2 * offset;

               // Here we find the values of the image at (i', j', k').
               // Interpolation is necessary because (i', j', k') are not
               // exactly integers.             
               const float g = linearInterpolator(
                  ip,
                  jp,
                  kp,
                  image,
                  dim.nx,
                  dim.ny,
                  dim.nz,
                  np
               );

               // If the reflected point also lies inside the object,
               // accumulate the statistics required for NCC.
               if(g != 0.0f)
               {
                  Sfg += f * g;
                  Sgg += g * g;
                  Sg += g;
                  Sff += f * f;
                  Sf += f;

                  numPairs++;
               }
            }

            q++;
         }
      }
   }

   if(numPairs == 0)
   {
      return 0.0f;
   }

   const double varF = Sff - Sf * Sf / numPairs;
   const double varG = Sgg - Sg * Sg / numPairs;

   if(varF <= DBL_EPSILON || varG <= DBL_EPSILON)
   {
      return 0.0f;
   }

   const double cov = Sfg - Sf * Sg / numPairs;

   float ncc = static_cast<float>(cov / std::sqrt(varF * varG));

   // Clamping the ncc because roundoff occasionally produces numbers
   // very slightly outside the range [-1, 1].
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
