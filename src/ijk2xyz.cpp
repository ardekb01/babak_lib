#include "babak_lib.h"

bool xyz2ijk(float *T,
             int nx,
             int ny,
             int nz,
             float dx,
             float dy,
             float dz)
{
   if (T == nullptr)
   {
      return false;
   }

   if (nx <= 0 || ny <= 0 || nz <= 0 ||
       dx <= 0.0f || dy <= 0.0f || dz <= 0.0f)
   {
      return false;
   }

   T[0]  = 1.0f / dx;
   T[1]  = 0.0f;
   T[2]  = 0.0f;
   T[3]  = (nx - 1.0f) / 2.0f;

   T[4]  = 0.0f;
   T[5]  = 1.0f / dy;
   T[6]  = 0.0f;
   T[7]  = (ny - 1.0f) / 2.0f;

   T[8]  = 0.0f;
   T[9]  = 0.0f;
   T[10] = 1.0f / dz;
   T[11] = (nz - 1.0f) / 2.0f;

   T[12] = 0.0f;
   T[13] = 0.0f;
   T[14] = 0.0f;
   T[15] = 1.0f;

   return true;
}

bool ijk2xyz(float *T,
             int nx,
             int ny,
             int nz,
             float dx,
             float dy,
             float dz)
{
   if (T == nullptr)
   {
      return false;
   }

   if (nx <= 0 || ny <= 0 || nz <= 0 ||
       dx <= 0.0f || dy <= 0.0f || dz <= 0.0f)
   {
      return false;
   }

   T[0]  = dx;
   T[1]  = 0.0f;
   T[2]  = 0.0f;
   T[3]  = -dx * (nx - 1.0f) / 2.0f;

   T[4]  = 0.0f;
   T[5]  = dy;
   T[6]  = 0.0f;
   T[7]  = -dy * (ny - 1.0f) / 2.0f;

   T[8]  = 0.0f;
   T[9]  = 0.0f;
   T[10] = dz;
   T[11] = -dz * (nz - 1.0f) / 2.0f;

   T[12] = 0.0f;
   T[13] = 0.0f;
   T[14] = 0.0f;
   T[15] = 1.0f;

   return true;
}

bool xyz2ijk(float *T, const DIM &dim)
{
   if (T == nullptr)
   {
      return false;
   }

   if (dim.dx <= 0.0f ||
       dim.dy <= 0.0f ||
       dim.dz <= 0.0f ||
       dim.nx <= 0 ||
       dim.ny <= 0 ||
       dim.nz <= 0)
   {
      return false;
   }

   T[0]  = 1.0f / dim.dx;
   T[1]  = 0.0f;
   T[2]  = 0.0f;
   T[3]  = (dim.nx - 1.0f) / 2.0f;

   T[4]  = 0.0f;
   T[5]  = 1.0f / dim.dy;
   T[6]  = 0.0f;
   T[7]  = (dim.ny - 1.0f) / 2.0f;

   T[8]  = 0.0f;
   T[9]  = 0.0f;
   T[10] = 1.0f / dim.dz;
   T[11] = (dim.nz - 1.0f) / 2.0f;

   T[12] = 0.0f;
   T[13] = 0.0f;
   T[14] = 0.0f;
   T[15] = 1.0f;

   return true;
}

bool ijk2xyz(float *T, const DIM &dim)
{
   if (T == nullptr)
   {
      return false;
   }

   if (dim.dx <= 0.0f ||
       dim.dy <= 0.0f ||
       dim.dz <= 0.0f ||
       dim.nx <= 0 ||
       dim.ny <= 0 ||
       dim.nz <= 0)
   {
      return false;
   }

   T[0]  = dim.dx;
   T[1]  = 0.0f;
   T[2]  = 0.0f;
   T[3]  = -dim.dx * (dim.nx - 1.0f) / 2.0f;

   T[4]  = 0.0f;
   T[5]  = dim.dy;
   T[6]  = 0.0f;
   T[7]  = -dim.dy * (dim.ny - 1.0f) / 2.0f;

   T[8]  = 0.0f;
   T[9]  = 0.0f;
   T[10] = dim.dz;
   T[11] = -dim.dz * (dim.nz - 1.0f) / 2.0f;

   T[12] = 0.0f;
   T[13] = 0.0f;
   T[14] = 0.0f;
   T[15] = 1.0f;

   return true;
}
