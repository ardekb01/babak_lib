#include "babak_lib.h"
#include "set_dim.h"
#include "ijk2xyz.h"

// Convert P from (x, y, z) to (i, j, k) coordinates.
// P is a 3 x n matrix. Each column is a point in the
// (x, y, z) coordinate system. Each point is converted
// in place to the (i, j, k) coordinate system.
// Core operation: P -> point -> matrix multiplication -> result -> P
bool convert_to_ijk(float *P, int n, SHORTIM im)
{
   // Validate P.
   if (P == nullptr)
   {
      return false;
   }

   // Validate n.
   if (n <= 0)
   {
      return false;
   }

   float T[16];
   float point[4];
   float result[4];
   DIM dim;

   set_dim(dim, im);

   if (!xyz2ijk(T, dim))
   {
      return false;
   }

   for (int i = 0; i < n; i++)
   {
      // Note that P is a 3 x n matrix.
      point[0] = P[i];
      point[1] = P[n + i];
      point[2] = P[2 * n + i];
      point[3] = 1.0f;

      multi(T, 4, 4, point, 4, 1, result);

      P[i]         = result[0];
      P[n + i]     = result[1];
      P[2 * n + i] = result[2];
   }

   return true;
}

// Convert P from (i, j, k) to (x, y, z) coordinates.
// P is a 3 x n matrix. Each column is a point in the
// (i, j, k) coordinate system. Each point is converted
// in-place to the (x, y, z) coordinate system.
// Core operation: P -> point -> matrix multiplication -> result -> P
bool convert_to_xyz(float *P, int n, SHORTIM im)
{
   // Validate P.
   if (P == nullptr)
   {
      return false;
   }

   // Validate n.
   if (n <= 0)
   {
      return false;
   }

   float T[16];
   float point[4];
   float result[4];
   DIM dim;

   set_dim(dim, im);

   if (!ijk2xyz(T, dim))
   {
      return false;
   }

   for (int i = 0; i < n; i++)
   {
      // Note that P is a 3 x n matrix.
      point[0] = P[i];
      point[1] = P[n + i];
      point[2] = P[2 * n + i];
      point[3] = 1.0f;

      multi(T, 4, 4, point, 4, 1, result);

      P[i] = result[0];
      P[n + i] = result[1];
      P[2 * n + i] = result[2];
   }

   return true;
}

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
