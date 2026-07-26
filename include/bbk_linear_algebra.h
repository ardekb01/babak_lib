#ifndef BBK_LINEAR_ALGEBRA_H 
#define BBK_LINEAR_ALGEBRA_H 

#define BBKTOL 1.0e-12

bool set_to_I(float *A, int n);
bool vectorNorm(const float *x, int n, double &norm);
bool normalizeVector(float *x, int n);
bool rotationMatrix(float *R, float alpha, float x, float y, float z);
bool rotationMatrix(float *R,
            float cosAlpha,
            float sinAlpha,
            float x,
            float y,
            float z);

float det3x3(const float *A);
double det3x3(const double *A);

float det4x4(const float *A);
double det4x4(const double *A);

float *inv2x2(const float *A);
double *inv2x2(const double *A);

float *inv3x3(const float *A);
double *inv3x3(const double *A);

// Transposes an (nRows × nCols) matrix stored in row-major order.
// The matrix is transposed in place using temporary storage.
// After the operation, the matrix has dimensions (nCols × nRows).
bool transpose_matrix(float *A, int nRows, int nCols);

// Transposes an (nRows × nCols) matrix stored in row-major order.
// The result is stored in transA and has dimensions (nCols × nRows).
bool transpose_matrix(const float *A,
                      int nRows,
                      int nCols,
                      float *transA);

// Returns a newly allocated transpose of an (nRows × nCols) matrix
// stored in row-major order. The returned matrix has dimensions
// (nCols × nRows). The caller is responsible for freeing the memory.
float *transpose(const float *A, int nRows, int nCols);
double *transpose(const double *A, int nRows, int nCols);

float *inv4(float *A);

#endif
