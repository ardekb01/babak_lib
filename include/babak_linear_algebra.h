#ifndef BABAK_LINEAR_ALGEBRA_H 
#define BABAK_LINEAR_ALGEBRA_H 

bool set_to_I(float *A, int n);
bool vectorNorm(const float *x, int n, double &norm);
bool normalizeVector(float *x, int n);
bool rotate(float *R, float alpha, float x, float y, float z);
bool rotate(float *R,
            float cosAlpha,
            float sinAlpha,
            float x,
            float y,
            float z);

#endif
