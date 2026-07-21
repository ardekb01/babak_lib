#include <math.h>

/*
   This function determines the general direction of a vector (x, y, z)
   defined in the RAS system.

   Possible directions are:
      Right, Left, Anterior, Posterior, Superior, or Inferior.

   Input:
      (x, y, z) - A vector defined in the RAS system.

   Output:
      One of six characters: {R, L, A, P, S, I}.

   If two or more components have equal maximum magnitudes,
   the priority is x, then y, then z.
*/

char directionCode(float x, float y, float z)
{
   float absX;
   float absY;
   float absZ;

   // Absolute values of x, y, and z.
   absX = fabsf(x);
   absY = fabsf(y);
   absZ = fabsf(z);

   if (absX >= absY && absX >= absZ)
   {
      if (x > 0.0f)
      {
         return 'R';
      }
      else
      {
         return 'L';
      }
   }
   else if (absY >= absX && absY >= absZ)
   {
      if (y > 0.0f)
      {
         return 'A';
      }
      else
      {
         return 'P';
      }
   }
   else
   {
      if (z > 0.0f)
      {
         return 'S';
      }
      else
      {
         return 'I';
      }
   }
}
