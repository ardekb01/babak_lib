#ifndef _linearInterpolator_h
#include <cstddef>

template<class TYPE> float linearInterpolator(float x, float y, TYPE *array, int nx, int ny)
{
   int i, j, n;
   float u, uu;
   float v1, v2;

   i = (int)(x);
   j = (int)(y);

   if (i < 0 || i > (nx - 2) || j < 0 || j > (ny - 2))
   {
      return(0.0);
   }

   u = x - i;
   if (u < 0.0)
      u = 0.0;

   uu = 1.0 - u;

   n = j * nx + i;
   v1 = array[n] * uu + array[n + 1] * u;
   v2 = array[n + nx] * uu + array[n + nx + 1] * u;

   u = y - j;
   if (u < 0.0)
      u = 0.0;

   uu = 1.0 - u;

   return(v1 * uu + v2 * u);
}

template<class TYPE> float linearInterpolator(float x, float y, float z,
                                              TYPE *array, int nx, int ny,
                                              int nz, int np)
{
   int i, j, k, n;
   float u, uu;
   float v1, v2, v3, v4;
   float w1, w2;

   i = (int)(x);
   j = (int)(y);
   k = (int)(z);

   if (i < 0 || i > (nx - 2) || j < 0 || j > (ny - 2))
   {
      return(0.0);
   }

   if (k >= 0 && k < (nz - 1))
   {
      u = x - i;
      if (u < 0.0)
         u = 0.0;

      uu = 1.0 - u;

      n = k * np + j * nx + i;

      v1 = array[n] * uu + array[n + 1] * u;
      v2 = array[n + nx] * uu + array[n + nx + 1] * u;
      v3 = array[n + np] * uu + array[n + np + 1] * u;
      v4 = array[n + np + nx] * uu + array[n + np + nx + 1] * u;

      u = y - j;
      if (u < 0.0)
         u = 0.0;

      uu = 1.0 - u;

      w1 = v1 * uu + v2 * u;
      w2 = v3 * uu + v4 * u;

      u = z - k;
      if (u < 0.0)
         u = 0.0;

      return(w1 * (1.0 - u) + w2 * u);
   }

   if (k == (nz - 1))
   {
      u = x - i;
      if (u < 0.0)
         u = 0.0;

      uu = 1.0 - u;

      n = k * np + j * nx + i;
      v1 = array[n] * uu + array[n + 1] * u;

      n = k * np + (j + 1) * nx + i;
      v2 = array[n] * uu + array[n + 1] * u;

      u = y - j;
      if (u < 0.0)
         u = 0.0;

      return(v1 * (1.0 - u) + v2 * u);
   }

   return(0.0);
}

template<class TYPE> float nearestNeighbor(float x, float y, float z,
                                           TYPE *array, int nx, int ny,
                                           int nz, int np)
{
   int i, j, k;

   i = (int)(x + 0.5);
   j = (int)(y + 0.5);
   k = (int)(z + 0.5);

   if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
   {
      return((float)(array[np * k + nx * j + i]));
   }
   else
   {
      return(0.0);
   }
}

#define _linearInterpolator_h
#endif
