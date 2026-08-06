#include <cmath>
#include <cfloat>
#include "babak_lib.h"

// have to limit the search considering the fact that the input image will be almost PIL
// The input image must be in PIL orientation.
bool findInitialNormalVector(short *image, const DIM &dim, float &A, float &B, float &C)
{
   if(image == nullptr)
      return false;

   if(dim.nx <= 0 || dim.ny <= 0 || dim.nz <= 0)
      return false;
  
   if(dim.dx <= 0.0f || dim.dy <= 0.0f || dim.dz <= 0.0f)
      return false;

   // compute x_cm,y_cm,z_cm the coordinates of the image "center of gravity"
   // in (mm) with respect to the FOV center as origin.
   float x_cm, y_cm, z_cm;
   if( !intensity_weighted_centroid(image, dim.nx, dim.ny, dim.nz,
      dim.dx, dim.dy, dim.dz,
      x_cm, y_cm, z_cm) )
   {
      return false;
   }

   // Ensure that the plane does not pass through the origin.
   if( std::abs(z_cm) <=  FLT_EPSILON)
   {
      z_cm = 0.1f;
   }

   float dum;
   float a, b, c;    // direction cosines  ax + by + cz = d 
   float d;

   // Ax + By + Cz = 1
   float A1, B1, C1;

   float cc;
   float ccmax = 0.0f;

   /* Coordinates of the input image "center of gravity" in mm.
      Origin is taken to be the center of the image volume. */

   double pi;
   double phi0 = 20.0;
   double delphi = 2.0;
   int nrings;                 // number of rings
   double ringlength;          // length of a ring
   double *cumulativelength;
   double *theta;
   double *phi;
   double totallength;         // sum of all ring lengths
   double arclength;
   int N;                      // number of samples

   pi = 4.0 * std::atan(1.0);

   // convert from degrees to radians
   phi0 = pi * phi0 / 180.0;
   delphi = pi * delphi / 180.0;

   nrings = static_cast<int>(std::ceil(phi0 / delphi)) + 1;

   cumulativelength = (double *)calloc(nrings, sizeof(double));
  
   if(cumulativelength == nullptr)
      return false;

   totallength = 0.0;
   for(int i = 0; i < nrings; i++)
   {
      ringlength = 2 * pi * sin(i * delphi);

      totallength += ringlength;

      cumulativelength[i] = totallength;
   }

   N = static_cast<int>(std::floor(totallength / delphi)) + 1;

   theta = (double *)calloc(N, sizeof(double));
   
   if(theta == nullptr)
   {
      free(cumulativelength);
      return false;
   }

   phi = (double *)calloc(N, sizeof(double));

   if(phi == nullptr)
   {
      free(cumulativelength);
      free(theta);
      return false;
   }

   theta[0] = 0.0;
   phi[0] = 0.0;

   for(int i = 1; i < N; i++)
   {
      arclength = i * delphi;

      for(int j = 1; j < nrings; j++)
      {
         if(arclength <= cumulativelength[j] &&
            arclength > cumulativelength[j - 1])
         {
            arclength -= cumulativelength[j - 1];
            phi[i] = j * delphi;
            theta[i] = arclength / sin(phi[i]);
            break;
         }
      }
   }

   //printf("\n******x_cm=%7.3f y_cm=%7.3f z_cm=%7.3f (mm)\n",x_cm,y_cm,z_cm);
   //printf("i = %f\n", (x_cm + dim.dx*(dim.nx-1.0)/2.0)/1.5 );
   //printf("j = %f\n", (y_cm + dim.dy*(dim.ny-1.0)/2.0)/0.859375 );
   //printf("k = %f\n", (z_cm + dim.dz*(dim.nz-1.0)/2.0)/0.859375 );

   A1 = 0.0f;
   B1 = 0.0f;
   C1 = 1.0f / z_cm;

   for(float z = z_cm - 5.0f; z <= z_cm + 5.0f; z += 1.0f)
   {
      for(int i = 0; i < N; i++)
      {
         /* The samples theta and phi define a direction in space. Find the
         unit vector (a,b,c) in that direction. */
         a = (float)(std::sin(phi[i]) * std::cos(theta[i]));
         b = (float)(std::sin(phi[i]) * std::sin(theta[i]));
         c = (float)std::cos(phi[i]);

         d = a * x_cm + b * y_cm + c * z;

         /* make sure d is non-negative */
         if(d < 0.0f)
         {
            a *= -1.0f;
            b *= -1.0f;
            c *= -1.0f;
            d *= -1.0f;
         }

         /* find the cross-correlation between image and its reflection
         about the plane ax+by+cz=d */

         if(d > FLT_EPSILON)
         {
            cc = symm_objective_func(image, dim, a / d, b / d, c / d);

            if(cc > ccmax)
            {
               ccmax = cc;
               A1 = a / d;
               B1 = b / d;
               C1 = c / d;
            }
         }
      }
   }

   dum = (float)std::sqrt((double)A1 * A1 + B1 * B1 + C1 * C1);
   a = A1 / dum;
   b = B1 / dum;
   c = C1 / dum;
   d = 1.0f / dum;

   //printf("\nInitial guess:");
   //printf("\nplane of symmetry: (%7.3f,%7.3f,%7.3f).(x,y,z) = %7.3f", a,b,c,d);
   //printf("\ncross correlation = %6.4f\n",ccmax);

   cc = optimizeNormalVector(image, dim, A1, B1, C1);

   //dum=(float)std::sqrt((double)A1*A1 + B1*B1 + C1*C1 );
   //a=A1/dum; b=B1/dum; c=C1/dum; d=1./dum;
   //printf("\nRefined initial guess:");
   //printf("\nplane of symmetry: (%7.3f,%7.3f,%7.3f).(x,y,z) = %7.3f", a,b,c,d);
   //printf("\ncross correlation = %6.4f\n",cc);

   A = A1;
   B = B1;
   C = C1;

   free(cumulativelength);
   free(theta);
   free(phi);

   return true;
}
