#include <cmath>
#include <cstdlib>
#include "bbk_linear_algebra.h"

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

// Uses Rodrigues' formula to produce a 4x4 transformation matrix for rotating
// a point by an angle alpha about the (x, y, z) axis.
bool rotationMatrix(float *R,
            float cosAlpha,
            float sinAlpha,
            float x,
            float y,
            float z)
{
   if (!set_to_I(R, 4))
   {
      return false;
   }

   if (x == 0.0f && y == 0.0f && z == 0.0f)
   {
      return false;
   }

   // Compute the Euclidean norm of the rotation axis.
   double norm = std::hypot(std::hypot((double)x, (double)y),
                            (double)z);

   // Find the unit vector (ax, ay, az) in the direction of (x, y, z).
   double ax = (double)x / norm;
   double ay = (double)y / norm;
   double az = (double)z / norm;

   R[0] = (float)(ax * ax + cosAlpha - cosAlpha * ax * ax);
   R[1] = (float)(ax * ay - cosAlpha * ax * ay - sinAlpha * az);
   R[2] = (float)(ax * az - cosAlpha * ax * az + sinAlpha * ay);

   R[4] = (float)(ay * ax - cosAlpha * ay * ax + sinAlpha * az);
   R[5] = (float)(ay * ay + cosAlpha - cosAlpha * ay * ay);
   R[6] = (float)(ay * az - cosAlpha * ay * az - sinAlpha * ax);

   R[8] = (float)(az * ax - cosAlpha * az * ax - sinAlpha * ay);
   R[9] = (float)(az * ay - cosAlpha * az * ay + sinAlpha * ax);
   R[10] = (float)(az * az + cosAlpha - cosAlpha * az * az);

   return true;
}

// Uses Rodrigues' formula to produce a 4x4 transformation matrix for rotating
// a point by an angle alpha about the (x, y, z) axis.
bool rotationMatrix(float *R, float alpha, float x, float y, float z)
{
   if (!set_to_I(R, 4))
      return false;

   if (x == 0.0f && y == 0.0f && z == 0.0f)
   {
      return false;
   }

   double cosAlpha = std::cos((double)alpha);
   double sinAlpha = std::sin((double)alpha);

   // Compute the Euclidean norm of the rotation axis.
   double norm = std::hypot(std::hypot((double)x, (double)y),
                            (double)z);

   // Find the unit vector (ax, ay, az) in the direction of (x, y, z).
   double ax = (double)x / norm;
   double ay = (double)y / norm;
   double az = (double)z / norm;

   R[0] = (float)(cosAlpha + ax * ax * (1.0 - cosAlpha));
   R[1] = (float)(ax * ay - cosAlpha * ax * ay - sinAlpha * az);
   R[2] = (float)(ax * az - cosAlpha * ax * az + sinAlpha * ay);

   R[4] = (float)(ay * ax - cosAlpha * ay * ax + sinAlpha * az);
   R[5] = (float)(cosAlpha + ay * ay * (1.0 - cosAlpha));
   R[6] = (float)(ay * az - cosAlpha * ay * az - sinAlpha * ax);

   R[8] = (float)(az * ax - cosAlpha * az * ax - sinAlpha * ay);
   R[9] = (float)(az * ay - cosAlpha * az * ay + sinAlpha * ax);
   R[10] = (float)(cosAlpha + az * az * (1.0 - cosAlpha));

   return true;
}

// Computes the determinant of a 3×3 matrix stored in row-major order.
float det3x3(const float *A)
{
   return A[0] * A[4] * A[8] +
          A[1] * A[5] * A[6] +
          A[2] * A[3] * A[7] -
          A[2] * A[4] * A[6] -
          A[0] * A[5] * A[7] -
          A[1] * A[3] * A[8];
}

double det3x3(const double *A)
{
   return A[0] * A[4] * A[8] +
          A[1] * A[5] * A[6] +
          A[2] * A[3] * A[7] -
          A[2] * A[4] * A[6] -
          A[0] * A[5] * A[7] -
          A[1] * A[3] * A[8];
}

// Computes the determinant of a 4×4 matrix stored in row-major order.
float det4x4(const float *A)
{
   float det = 0.0f;
   float B[9];

   B[0] = A[5];
   B[1] = A[6];
   B[2] = A[7];
   B[3] = A[9];
   B[4] = A[10];
   B[5] = A[11];
   B[6] = A[13];
   B[7] = A[14];
   B[8] = A[15];

   det += A[0] * det3x3(B);

   B[0] = A[4];
   B[1] = A[6];
   B[2] = A[7];
   B[3] = A[8];
   B[4] = A[10];
   B[5] = A[11];
   B[6] = A[12];
   B[7] = A[14];
   B[8] = A[15];

   det -= A[1] * det3x3(B);

   B[0] = A[4];
   B[1] = A[5];
   B[2] = A[7];
   B[3] = A[8];
   B[4] = A[9];
   B[5] = A[11];
   B[6] = A[12];
   B[7] = A[13];
   B[8] = A[15];

   det += A[2] * det3x3(B);

   B[0] = A[4];
   B[1] = A[5];
   B[2] = A[6];
   B[3] = A[8];
   B[4] = A[9];
   B[5] = A[10];
   B[6] = A[12];
   B[7] = A[13];
   B[8] = A[14];

   det -= A[3] * det3x3(B);

   return det;
}

// Computes the determinant of a 4×4 matrix stored in row-major order.
double det4x4(const double *A)
{
   double det = 0.0;
   double B[9];

   B[0] = A[5];
   B[1] = A[6];
   B[2] = A[7];
   B[3] = A[9];
   B[4] = A[10];
   B[5] = A[11];
   B[6] = A[13];
   B[7] = A[14];
   B[8] = A[15];

   det += A[0] * det3x3(B);

   B[0] = A[4];
   B[1] = A[6];
   B[2] = A[7];
   B[3] = A[8];
   B[4] = A[10];
   B[5] = A[11];
   B[6] = A[12];
   B[7] = A[14];
   B[8] = A[15];

   det -= A[1] * det3x3(B);

   B[0] = A[4];
   B[1] = A[5];
   B[2] = A[7];
   B[3] = A[8];
   B[4] = A[9];
   B[5] = A[11];
   B[6] = A[12];
   B[7] = A[13];
   B[8] = A[15];

   det += A[2] * det3x3(B);

   B[0] = A[4];
   B[1] = A[5];
   B[2] = A[6];
   B[3] = A[8];
   B[4] = A[9];
   B[5] = A[10];
   B[6] = A[12];
   B[7] = A[13];
   B[8] = A[14];

   det -= A[3] * det3x3(B);

   return det;
}

// Computes the inverse of a 2×2 matrix stored in row-major order.
// Returns a newly allocated matrix, or nullptr if A is null, singular,
// or memory allocation fails.
// Caller is responsible for freeing the memory in the returned matrix
// using free().
float *inv2x2(const float *A)
{
   if (A == nullptr)
   {
      return nullptr;
   }

   float detA = A[0] * A[3] -
                A[1] * A[2];

   if (fabsf(detA) < BBKTOL)
   {
      return nullptr;
   }

   float *invA = (float *)calloc(4, sizeof(float));

   if (invA == nullptr)
   {
      return nullptr;
   }

   invA[0] = A[3] / detA;
   invA[1] = -A[1] / detA;
   invA[2] = -A[2] / detA;
   invA[3] = A[0] / detA;

   return invA;
}

double *inv2x2(const double *A)
{
   if (A == nullptr)
   {
      return nullptr;
   }

   double detA = A[0] * A[3] -
                A[1] * A[2];

   if (fabs(detA) < BBKTOL)
   {
      return nullptr;
   }

   double *invA = (double *)calloc(4, sizeof(double));

   if (invA == nullptr)
   {
      return nullptr;
   }

   invA[0] = A[3] / detA;
   invA[1] = -A[1] / detA;
   invA[2] = -A[2] / detA;
   invA[3] = A[0] / detA;

   return invA;
}
