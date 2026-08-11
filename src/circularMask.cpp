#include <cstddef>
#include <cstdlib>
#include <cmath>
#include "bbk.h"

// center[0] = physical x-coordinate relative to image center
// center[1] = physical y-coordinate relative to image center
char *circularMask(short *mask, 
                   const DIM &dim, 
                   const float *center, 
                   const double radius)
{
   if(dim.nx <= 0 || dim.ny <= 0 || dim.nz <= 0)
      return nullptr;

   if(dim.dx <= 0.0f || dim.dy <= 0.0f)
      return nullptr;

   if(radius < 0.0)
      return nullptr;

   if(mask == nullptr || center == nullptr)
      return nullptr;

   const size_t np = static_cast<size_t>(dim.nx) *
                     static_cast<size_t>(dim.ny);

   char *circle = static_cast<char *>(malloc(np));

   if( circle == nullptr )
   {
      return nullptr;
   }

   const double xc = (dim.nx-1.0)/2.0;
   const double yc = (dim.ny-1.0)/2.0;

   for(int j=0; j<dim.ny; j++)
   {
      const double r2 = (j - yc)*dim.dy - center[1];

      for(int i=0; i<dim.nx; i++)
      {
         const double r1 = (i - xc)*dim.dx - center[0];

         const size_t index = static_cast<size_t>(j) *
                              static_cast<size_t>(dim.nx) +
                              static_cast<size_t>(i);

         const double r  = std::hypot(r1 , r2);

         if(r >= radius) 
         {
            circle[index]=0;

            for(int k=0; k<dim.nz; k++)
            {
               mask[static_cast<size_t>(k) * np + index] = 0;
            }
         }
         else
         {	
            circle[index]=1;
         }
      }
   }

   return(circle);
}
