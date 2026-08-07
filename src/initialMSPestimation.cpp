#include <cmath>
#include <vector>
#include <cfloat>
#include "babak_lib.h"

// have to limit the search considering the fact that the input image will be almost PIL
// The input image must be in PIL orientation.
// Plane parameters are represented as Ax + By + Cz = 1.
bool initialMSPestimation(const short *image, const DIM &dim, float &A, float &B, float &C)
{
   if(image == nullptr || 
     dim.nx <= 0 || dim.ny <= 0 || dim.nz <= 0 ||
     dim.dx <= 0.0f || dim.dy <= 0.0f || dim.dz <= 0.0f)
   {
      return false;
   }

   // Compute x_cm,y_cm,z_cm the coordinates of the intensity-weighted image 
   // centroid in (mm) with respect to the FOV center as origin.
   float x_cm, y_cm, z_cm;
   if( !intensity_weighted_centroid(image, dim.nx, dim.ny, dim.nz,
      dim.dx, dim.dy, dim.dz,
      x_cm, y_cm, z_cm) )
   {
      return false;
   }

   // Maximum angular deviation from the superior-inferior axis.
   // Assumes the input image is approximately PIL-oriented.
   constexpr double kMaxPhiDeg = 20.0;

   constexpr double kPhiStepDeg = 2.0;
   constexpr float kMinimumOffset = 0.1f;
   constexpr double kPi = 3.14159265358979323846;
   constexpr double kDegToRad = kPi / 180.0;
   constexpr int kHalfSearchRange = 5;

   // Ensure that the plane does not pass through the origin.
   if( std::abs(z_cm) <=  FLT_EPSILON)
   {
      z_cm = kMinimumOffset;
   }

   float a, b, c;    // direction cosines  ax + by + cz = d 
   float d;
   float ccmax = -1.0f;
   double phi0 = kMaxPhiDeg;
   double delphi = kPhiStepDeg;
   int nrings;                 // number of rings
   double ringLength;          // length of a ring
   double totalLength;         // sum of all ring lengths
   double arclength;
   int N;                      // number of samples

   // Convert from degrees to radians.
   phi0 *= kDegToRad;
   delphi *= kDegToRad;

   nrings = static_cast<int>(std::ceil(phi0 / delphi)) + 1;

   std::vector<double> cumulativeLength(nrings);

   totalLength = 0.0;
   for(int i = 0; i < nrings; i++)
   {
      ringLength = 2 * kPi * std::sin(i * delphi);

      totalLength += ringLength;

      cumulativeLength[i] = totalLength;
   }

   N = static_cast<int>(std::floor(totalLength / delphi)) + 1;

   std::vector<double> theta(N);
   std::vector<double> phi(N);

   theta[0] = 0.0;
   phi[0] = 0.0;

   for(int i = 1; i < N; i++)
   {
      arclength = i * delphi;

      for(int j = 1; j < nrings; j++)
      {
         if(arclength <= cumulativeLength[j] &&
            arclength > cumulativeLength[j - 1])
         {
            arclength -= cumulativeLength[j - 1];
            phi[i] = j * delphi;
            theta[i] = arclength / std::sin(phi[i]);
            break;
         }
      }
   }

   std::vector<float> dirX(N);
   std::vector<float> dirY(N);
   std::vector<float> dirZ(N);

   for(int i = 0; i < N; i++)
   {
      const float sinPhi   = static_cast<float>(std::sin(phi[i]));
      const float cosPhi   = static_cast<float>(std::cos(phi[i]));
      const float sinTheta = static_cast<float>(std::sin(theta[i]));
      const float cosTheta = static_cast<float>(std::cos(theta[i]));

      dirX[i] = sinPhi * cosTheta;
      dirY[i] = sinPhi * sinTheta;
      dirZ[i] = cosPhi;
   }

   A = 0.0f;
   B = 0.0f;
   C = 1.0f / z_cm;

   for(int offset = -kHalfSearchRange; offset <= kHalfSearchRange; offset++)
   {
      const float z = z_cm + static_cast<float>(offset);

      for(int i = 0; i < N; i++)
      {
         // The samples theta and phi define a direction in space. Find the
         // unit vector (a,b,c) in that direction.
         a = dirX[i];
         b = dirY[i];
         c = dirZ[i];

         d = a * x_cm + b * y_cm + c * z;

         // make sure d is non-negative 
         if(d < 0.0f)
         {
            a *= -1.0f;
            b *= -1.0f;
            c *= -1.0f;
            d *= -1.0f;
         }

         // find the cross-correlation between image and its reflection
         // about the plane ax+by+cz=d 

         if(d > FLT_EPSILON)
         {
            const float invD = 1.0f / d;

            const float cc =
               symm_objective_func(image, dim,
                                   a * invD,
                                   b * invD,
                                   c * invD);

            if(cc > ccmax)
            {
               ccmax = cc;
               A = a * invD;
               B = b * invD;
               C = c * invD;
            }
         }
      }
   }

   if(ccmax < 0.0f)
      return false;

   return true;
}
