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

#endif
