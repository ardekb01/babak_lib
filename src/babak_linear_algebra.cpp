#include <math.h>
#include "babak_linear_algebra.h"

// Sets the n x n matrix A equal to the identity matrix.
bool set_to_I(float *A, int n)
{
   if (A == nullptr || n <= 0)
   {
      return false;
   }

   // Set all elements to zero.
   for (int i = 0; i < n * n; i++)
   {
      A[i] = 0.0f;
   }

   // Set the diagonal elements to one.
   for (int i = 0; i < n; i++)
   {
      A[n * i + i] = 1.0f;
   }

   return true;
}

// Normalize x in-place, that is: ||x||=1.0
bool normalizeVector(float *x, int n)
{
   double norm = 0.0;

   if (x == nullptr)
   {
      return false;
   }

   if (n <= 0)
   {
      return false;
   }

   if (!vectorNorm(x, n, norm))
   {
      return false;
   }

   if (norm == 0.0)
   {
      return false;
   }

   for (int i = 0; i < n; i++)
   {
      x[i] = static_cast<float>(x[i] / norm);
   }

   return true;
}

bool vectorNorm(const float *x, int n, double &norm)
{
   norm = 0.0;

   if (x == nullptr || n <= 0)
   {
      return false;
   }

   for (int i = 0; i < n; i++)
   {
      norm += x[i] * x[i];
   }

   norm = sqrt(norm);

   return true;
}
