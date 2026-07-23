#ifndef BBK_LINEAR_ALGEBRA_H 
#define BBK_LINEAR_ALGEBRA_H 

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

#endif
