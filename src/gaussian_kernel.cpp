/**************************************************
   program: gaussian_kernel.c
   Copyright (c) 1997 Babak A. Ardekani
   ALL RIGHTS RESERVED

   created: 25 June 1997
   revision: 26 June 1997 (BAA)
      gaussian_kernel() modified to handle the case where the input sd<=0
   revision: 3 August 2010 (BAA)
      removed the third input parameter (v) for printing
   revision: 4 July 2026 (BAA)
      modernized implementation with improved numerical robustness,
      overflow checking, and documentation (with assistance from ChatGPT).
      (BAA from 1997: "Assistance from who????")
***************************************************/

#include "gaussian_kernel.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

// factor 2.57 was chosen to make the area under the Gaussian approx. 0.99
// factor 1.65 would give an area of approx 0.90
#define GAUSSIAN_TRUNCATION_SIGMA 2.57

/**************************************************
Inputs 
   sd - standard deviation of the Gaussian IN UNITS OF PIXELS

Outputs
   n - dimension of the returned array
   returned pointer - array of dimension n (h[0],h[1],...,h[n-1])

Notes
   This function computes and returns an array of values proportional to:
   exp( -0.5*x^2/sd^2 ) ( x=0, 1, 2, ..., n-1 ).

   Since Gaussian curves are symmetric, only the right-half of the curve
   (including 0) is computed and returned (n points). The complete curve
   would have 2*n-1 points.

   The returned array is scaled such that the sum of all the coefficients
   (including the left-half of the curve) would be 1.

   The dimension n is chosen such that the area under the Gaussian curve is
   at least 0.99.

   The caller is responsible for freeing the returned array in 'h' using free().
***************************************************/

float *gaussian_kernel(const float sd, int *n)
{
   float *h;  /* the array to be returned */
   double var; /* constant variance */
   double nf;  /* normalization factor */

   if (n == NULL)
      return NULL;

   // Non-positive standard deviations produce a unit impulse kernel.
   if (sd <= 0.0f)
   {
      *n = 1;

      h = (float *)malloc(sizeof(float));

      if (h == NULL)
      {
         *n = 0;
         return NULL;
      }

      h[0] = 1.0f;

      return h;
   }

   var = (double)sd * sd;

   double size = ceil((double)sd * GAUSSIAN_TRUNCATION_SIGMA);

   if (size < 1.0)
   {
      size = 1.0;
   }

   // avoids integer overflow if sd is very large
   if (size > INT_MAX)
   {
      *n = 0;
      return NULL;
   }

   *n = (int)size;

   h = (float *)malloc(*n * sizeof(float));

   if (h == NULL)
   {
      *n = 0;
      return NULL;
   }

   nf = 0.0;

   for (int x = 0; x < (*n); x++)
   {
      h[x] = (float)exp(-0.5 * x * x / var);
      nf += h[x];
   }

   /* add the values of the left-half of the Gaussian curve to nf */
   nf = 2.0 * nf - h[0];

   /* normalize the Gaussian curve to 1 */
   for (int x = 0; x < (*n); x++)
   {
      h[x] = (float)(h[x]/nf);
   }

   return h;
}
