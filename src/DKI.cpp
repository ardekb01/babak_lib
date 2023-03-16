#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <babak_lib.h>

#define _DKI

////////////////////////////////////////////////////////////////////////////////
/////////////////////////Function prototypes////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void formA7(float *A, int n, float *u1, float *u2, float *u3, float *b, float *w);

float quadratic_form_1(float x0x0, float x0x1, float x0x2, float x1x1, float x1x2, float x2x2, float *v);
float quadratic_form_1(float x0, float x1, float x2, float *v);
float quadratic_form_1(float *x, float *v);
void vector_to_symmetric_matrix_form(float *u, float *D);
void symmetric_matrix_to_vector_form(float *D, float *u);
void cholesky_composition_1(float *v, float *u);
void cholesky_composition_2(float *L, float *D);
int cholesky_decomposition_1(float *u, float *v);
int cholesky_decomposition_2(float *D, float *L);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*******************Function prototypes**************************************/
void restricted_tensor_model(float *x, float *y);
void second_and_fourth_order_moments(float *x, float *y, float *M2, float *M4, int v);
void positive_definite_tensor_model(float *d1, float *d2, float *y);
/****************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/////////////////////////Function definitions///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// w: (nx1) weight vectors
void formA7(float *A, int n, float *u1, float *u2, float *u3, float *b, float *w)
{
   int row_offset;

   // i indexes rows of (n x 7) matrix A
   for(int i=0; i<n; i++)
   {
      row_offset = i*7; // skips to the begining of the i_th row (i=0,1,...,n-1)

      // set the 7 column enteries of A for the i_th row
      A[row_offset]   = 1.0*w[i]; 
      A[row_offset+1] = -b[i]*u1[i]*u1[i]*w[i];      // D11
      A[row_offset+2] = -2.0*b[i]*u1[i]*u2[i]*w[i];  // D12
      A[row_offset+3] = -2.0*b[i]*u1[i]*u3[i]*w[i];  // D13
      A[row_offset+4] = -b[i]*u2[i]*u2[i]*w[i];      // D22
      A[row_offset+5] = -2.0*b[i]*u2[i]*u3[i]*w[i];  // D23
      A[row_offset+6] = -b[i]*u3[i]*u3[i]*w[i];      // D33
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////

float estimate_K(float *A, float *x, float *y, int n)
{
   float Knum=0.0;
   float Kden=0.0;
   float *Ax;
   float K;
   float Ax2;

   Ax = (float *)calloc(n, sizeof(float));

   multi(A, n, 6 , x, 6, 1, Ax);

   for(int i=0; i<n; i++) 
   {
      Ax2 = (Ax[i]*Ax[i]);
      Knum += (y[i] - Ax[i])*Ax2;
      Kden += Ax2*Ax2;
   }
   K = Knum/Kden;

   free(Ax);

   return(K);
}

////////////////////////////////////////////////////////////////////////////////

void update_x(float *A, float *x, float *y, float K, int n)
{
   float *Ax;
   float e, ep, Aij;
   float Jp, Jpp;

   Ax = (float *)calloc(n, sizeof(float));

   for(int j=0; j<6; j++)
   {
      multi(A, n, 6 , x, 6, 1, Ax);

      Jp=Jpp=0.0;
      for(int i=0; i<n; i++)
      {
         e = y[i] - Ax[i] - K*Ax[i]*Ax[i];

         Aij = A[i*6+j];
         ep = -Aij*(1.0 + 2*K*Ax[i]);

         Jp += e*ep;
         Jpp += (-2.0*e*ep*K*Aij*Aij + ep*ep);
      }

      x[j] -= Jp/Jpp;
   }

   free(Ax);
}

////////////////////////////////////////////////////////////////////////////////

void update_v(float *u0, float *u1, float *u2, float *b, float *y, int n, float *v1, float *v2, float *v3, float scale)
{
   float Q1, Q2, Q3;
   float d1[6], d2[6], d3[6];
   float dJ1[6], dJ2[6], dJ3[6];
   float q0, q1, q2;
   float e;

   for(int i=0; i<6; i++) dJ1[i]=dJ2[i]=dJ3[i]=0.0;

   cholesky_composition_1(v1, d1);
   cholesky_composition_1(v2, d2);
   cholesky_composition_1(v3, d3);

   for(int i=0; i<n; i++)
   {
      q0 = u0[i]*sqrtf(b[i]/scale);
      q1 = u1[i]*sqrtf(b[i]/scale);
      q2 = u2[i]*sqrtf(b[i]/scale);

      Q1 = quadratic_form_1(q0, q1, q2, d1);
      Q2 = quadratic_form_1(q0, q1, q2, d2);
      Q3 = quadratic_form_1(q0, q1, q2, d3);

      e = y[i] + Q1 - Q2*Q3;

      dJ1[0] += e*( 2.0*(v1[0]*q0*q0 + v1[1]*q0*q1 + v1[2]*q0*q2) );
      dJ1[1] += e*( 2.0*(v1[0]*q0*q1 + v1[1]*q1*q1 + v1[2]*q1*q2) );
      dJ1[2] += e*( 2.0*(v1[0]*q0*q2 + v1[1]*q1*q2 + v1[2]*q2*q2) );
      dJ1[3] += e*( 2.0*(v1[3]*q1*q1 + v1[4]*q1*q2) );
      dJ1[4] += e*( 2.0*(v1[3]*q1*q2 + v1[4]*q2*q2) );
      dJ1[5] += e*( 2.0*v1[5]*q2*q2 );

      // dJ2[0] += e*( 2.0*(v2[0]*q0*q0 + v2[1]*q0*q1 + v2[2]*q0*q2) )*(-Q3);
      dJ2[1] += e*( 2.0*(v2[0]*q0*q1 + v2[1]*q1*q1 + v2[2]*q1*q2) )*(-Q3);
      dJ2[2] += e*( 2.0*(v2[0]*q0*q2 + v2[1]*q1*q2 + v2[2]*q2*q2) )*(-Q3);
      // dJ2[3] += e*( 2.0*(v2[3]*q1*q1 + v2[4]*q1*q2) )*(-Q3);
      dJ2[4] += e*( 2.0*(v2[3]*q1*q2 + v2[4]*q2*q2) )*(-Q3);
      // dJ2[5] += e*( 2.0*v2[5]*q2*q2 )*(-Q3);

      dJ3[0] += e*( 2.0*(v3[0]*q0*q0 + v3[1]*q0*q1 + v3[2]*q0*q2) )*(-Q2);
      dJ3[1] += e*( 2.0*(v3[0]*q0*q1 + v3[1]*q1*q1 + v3[2]*q1*q2) )*(-Q2);
      dJ3[2] += e*( 2.0*(v3[0]*q0*q2 + v3[1]*q1*q2 + v3[2]*q2*q2) )*(-Q2);
      dJ3[3] += e*( 2.0*(v3[3]*q1*q1 + v3[4]*q1*q2) )*(-Q2);
      dJ3[4] += e*( 2.0*(v3[3]*q1*q2 + v3[4]*q2*q2) )*(-Q2);
      dJ3[5] += e*( 2.0*v3[5]*q2*q2 )*(-Q2);
   }

   for(int i=0; i<6; i++)
   {
      v1[i] -= 0.01*dJ1[i];
      v2[i] -= 0.01*dJ2[i];
      v3[i] -= 0.01*dJ3[i];
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////

void update_v(float *u0, float *u1, float *u2, float *b, float *y, int n, float *v, float scale)
{
   float Q;
   float u[6];
   float dJ[6];
   float dQ[6];
   float q0, q1, q2;

   for(int i=0; i<6; i++) dJ[i]=0.0;

   cholesky_composition_1(v, u);

   for(int i=0; i<n; i++)
   {
      q0 = u0[i]*sqrtf(b[i]/scale);
      q1 = u1[i]*sqrtf(b[i]/scale);
      q2 = u2[i]*sqrtf(b[i]/scale);

      Q = y[i] + quadratic_form_1(q0, q1, q2, u);

      dQ[0] = 2.0*(v[0]*q0*q0 + v[1]*q0*q1 + v[2]*q0*q2);
      dQ[1] = 2.0*(v[0]*q0*q1 + v[1]*q1*q1 + v[2]*q1*q2);
      dQ[2] = 2.0*(v[0]*q0*q2 + v[1]*q1*q2 + v[2]*q2*q2);
      dQ[3] = 2.0*(v[3]*q1*q1 + v[4]*q1*q2);
      dQ[4] = 2.0*(v[3]*q1*q2 + v[4]*q2*q2);
      dQ[5] = 2.0*v[5]*q2*q2;

      for(int j=0; j<6; j++) dJ[j] += Q*dQ[j];
   }

   v[0] -= 0.01*dJ[0];
   v[1] -= 0.01*dJ[1];
   v[2] -= 0.01*dJ[2];
   v[3] -= 0.01*dJ[3];
   v[4] -= 0.01*dJ[4];
   v[5] -= 0.01*dJ[5];

   return;
}
////////////////////////////////////////////////////////////////////////////////

// Input 1: 3x1 vector x
// Input 2: 6x1 vector v representing a 3x3 symmetric matrix A such that:

//     v0 v1 v2
// A = v1 v3 v4
//     v2 v4 v5

// Ouput: returns x^T * A * x

float quadratic_form_1(float x0, float x1, float x2, float *v)
{
   float Q;

   Q = x0*x0*v[0] + 2.0*x0*x1*v[1] + 2.0*x0*x2*v[2] +
   x1*x1*v[3] + 2.0*x1*x2*v[4] + x2*x2*v[5];

   return(Q);
}

float quadratic_form_1(float x0x0, float x0x1, float x0x2, float x1x1, float x1x2, float x2x2, float *v)
{
   float Q;

   Q = x0x0*v[0] + 2.0*x0x1*v[1] + 2.0*x0x2*v[2] +
   x1x1*v[3] + 2.0*x1x2*v[4] + x2x2*v[5];

   return(Q);
}


////////////////////////////////////////////////////////////////////////////////

// Input 1: 3x1 vector x
// Input 2: 6x1 vector v representing a 3x3 symmetric matrix A such that:

//     v0 v1 v2
// A = v1 v3 v4
//     v2 v4 v5

// Ouput: returns x^T * A * x

float quadratic_form_1(float *x, float *v)
{
   float Q;

   Q = x[0]*x[0]*v[0] + 2.0*x[0]*x[1]*v[1] + 2.0*x[0]*x[2]*v[2] +
   x[1]*x[1]*v[3] + 2.0*x[1]*x[2]*v[4] + x[2]*x[2]*v[5];

   return(Q);
}

////////////////////////////////////////////////////////////////////////////////

// Input: 6x1 vector u
// Ouput: 3x3 symmetric matrix D such that:

//      u0  u1  u2
// D =  u1  u3  u4
//      u2  u4  u5

// Assumes a particular correspondence between elements of x and D
void vector_to_symmetric_matrix_form(float *u, float *D)
{
   D[0] = u[0];
   D[1] = D[3] = u[1];
   D[2] = D[6] = u[2];
   D[4] = u[3];
   D[5] = D[7] = u[4];
   D[8] = u[5];

   return;
}

// Assumes a particular correspondence between elements of x and D
void symmetric_matrix_to_vector_form(float *D, float *u)
{
   u[0] = D[0];
   u[1] = D[1];
   u[2] = D[2];
   u[3] = D[4];
   u[4] = D[5];
   u[5] = D[8];

   return;
}

////////////////////////////////////////////////////////////////////////////////

// Input: 6x1 vector v representing a 3x3 lower triangular matrix L

//     v0   0   0 
// L = v1  v3   0  
//     v2  v4  v5

// Outut: 6x1 vector u representing a 3x3 symmetric matrix D = L * L^T

//      u0  u1  u2
// D =  u1  u3  u4
//      u2  u4  u5

// Documented in Technical Notes

void cholesky_composition_1(float *v, float *u)
{
   u[0] = v[0]*v[0];
   u[1] = v[0]*v[1];
   u[2] = v[0]*v[2];
   u[3] = v[1]*v[1] + v[3]*v[3];
   u[4] = v[1]*v[2] + v[3]*v[4];
   u[5] = v[2]*v[2] + v[4]*v[4] + v[5]*v[5];

   return;
}

////////////////////////////////////////////////////////////////////////////////

// Input: 3x3 lower triangular matrix L

//     L0  0  0 
// L = L3 L4  0
//     L6 L7 L8

// Outut: 3x3 symmetric matrix D = L * L^T

// Documented in Technical Notes

void cholesky_composition_2(float *L, float *D)
{
   D[0] = L[0]*L[0];
   D[1] = L[0]*L[3];
   D[2] = L[0]*L[6];
   D[3] = D[1];
   D[4] = L[3]*L[3] + L[4]*L[4];
   D[5] = L[3]*L[6] + L[4]*L[7];
   D[6] = D[2];
   D[7] = D[5];
   D[8] = L[6]*L[6] + L[7]*L[7] + L[8]*L[8];

   return;
}

////////////////////////////////////////////////////////////////////////////////

// Input: 6x1 vector u representing 3x3 symmetric matrix D 

//      D0  D1  D2    u0 u1 u2
// D =  D3  D4  D5 =  u1 u3 u4
//      D6  D7  D8    u2 u4 u5

// Output: 6x1 vector v represnting 3x3 lower triangular matrix L

//      L0  0   0    v0 0  0
// L =  L3  L4  0  = v1 v2 0
//      L6  L7  L8   v3 v4 v5

// such that D = L * L^T

// returns 1 on error, zero otherwise
// Documented in Technical Notes

int cholesky_decomposition_1(float *u, float *v)
{
   float x;

   x=u[0];
   if(x<=0.0) return(1);
   v[0] = sqrtf(x);  

   v[1] = u[1]/v[0];    
   v[2] = u[2]/v[0];    

   x = u[3] - v[1]*v[1];
   if(x<=0.0) return(1);
   v[3] = sqrtf(x);   

   v[4] = (u[4] - v[1]*v[2])/v[3];

   x = u[5] - v[2]*v[2] - v[4]*v[4];
   if(x<=0.0) return(1);
   v[5] = sqrtf(x);

   return(0);
}

////////////////////////////////////////////////////////////////////////////////

// Input: 3x3 symmetric matrix D 

//      D0  D1  D2
// D =  D3  D4  D5
//      D6  D7  D8

// where D1=D3, D2=D6, and D5=D7.

// Output: 3x3 lower triangular matrix L

//      L0  L1  L2
// L =  L3  L4  L5
//      L6  L7  L8

// where L1=L2=L5=0, that is:

//      L0   0   0
// L =  L3  L4   0
//      L6  L7  L8

// such that D = L * L^T

//      L0*L0     L0*L3           L0*L6
// D =  L0*L3  L3*L3+L4*L4     L3*L6+L4*L7
//      L0*L6  L3*L6+L4*L7  L6*L6+L7*L7+L8*L8

// returns 1 on error, zero otherwise
// Documented in Technical Notes

int cholesky_decomposition_2(float *D, float *L)
{
   float x;

   x=D[0];
   if(x<=0.0) return(1);
   L[0] = sqrtf(x);  

   L[1] = 0.0;
   L[2] = 0.0;
   L[3] = D[1]/L[0];    

   x = D[4] - L[3]*L[3];
   if(x<=0.0) return(1);
   L[4] = sqrtf(x);   

   L[5] = 0.0;
   L[6] = D[2]/L[0];    
   L[7] = (D[5] - L[3]*L[6])/L[4];

   x = D[8] - L[6]*L[6] - L[7]*L[7];
   if(x<=0.0) return(1);
   L[8] = sqrtf(x);

   return(0);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// This function writes the elements of a general 4th order tensor (y) in terms 
// of the elements of the elements of a `restricted' 4th order tensor (x).
// See `Technical Notes' for details.

// x[0] = D11
// x[1] = D22
// x[2] = D33
// x[3] = D12
// x[4] = D13
// x[5] = D23

// y[0]  = D1111
// y[1]  = D2222
// y[2]  = D3333

// y[3]  = D1112
// y[4]  = D1113
// y[5]  = D2221
// y[6]  = D2223
// y[7]  = D3331
// y[8]  = D3332

// y[9]  = D1122
// y[10] = D1133
// y[11] = D2233

// y[12] = D1123
// y[13] = D2213
// y[14] = D3312

void restricted_tensor_model(float *x, float *y)
{
   y[0] = x[0]*x[0];
   y[1] = x[1]*x[1];
   y[2] = x[2]*x[2];

   y[3] = x[0]*x[3];
   y[4] = x[0]*x[4];
   y[5] = x[1]*x[3];
   y[6] = x[1]*x[5];
   y[7] = x[2]*x[4];
   y[8] = x[2]*x[5];

   y[9]  = (x[0]*x[1] + 2.0*x[3]*x[3])/3.0;
   y[10] = (x[0]*x[2] + 2.0*x[4]*x[4])/3.0;
   y[11] = (x[1]*x[2] + 2.0*x[5]*x[5])/3.0;

   y[12] = (x[0]*x[5] + 2.0*x[3]*x[4])/3.0;
   y[13] = (x[1]*x[4] + 2.0*x[3]*x[5])/3.0;
   y[14] = (x[2]*x[3] + 2.0*x[4]*x[5])/3.0;

   return;
}

// computes coefficients of (q^T D1 q)*(q^T D2 q)
// d1 and d2 are vector representations of D1 and D2
void positive_definite_tensor_model(float *d1, float *d2, float *y)
{
   y[0] = d1[0]*d2[0];
   y[1] = d1[3]*d2[3];
   y[2] = d1[5]*d2[5];

   y[3] = (d1[0]*d2[1] + d1[1]*d2[0])/2.0;
   y[4] = (d1[0]*d2[2] + d1[2]*d2[0])/2.0;
   y[5] = (d1[3]*d2[1] + d1[1]*d2[3])/2.0;
   y[6] = (d1[3]*d2[4] + d1[4]*d2[3])/2.0;
   y[7] = (d1[5]*d2[2] + d1[2]*d2[5])/2.0;
   y[8] = (d1[5]*d2[4] + d1[4]*d2[5])/2.0;

   y[9]  = (d1[0]*d2[3] + d1[3]*d2[0] + 4.0*d1[1]*d2[1])/6.0;
   y[10] = (d1[0]*d2[5] + d1[5]*d2[0] + 4.0*d1[2]*d2[2])/6.0;
   y[11] = (d1[3]*d2[5] + d1[5]*d2[3] + 4.0*d1[4]*d2[4])/6.0;

   y[12] = (d1[0]*d2[4] + d1[4]*d2[0] + 2.0*d1[1]*d2[2] + 2.0*d1[2]*d2[1])/3.0;
   y[13] = (d1[3]*d2[2] + d1[2]*d2[3] + 2.0*d1[1]*d2[4] + 2.0*d1[4]*d2[1])/3.0;
   y[14] = (d1[5]*d2[1] + d1[1]*d2[5] + 2.0*d1[2]*d2[4] + 2.0*d1[4]*d2[2])/3.0;

   return;
}

// M2[0] = M11 = 2*D11
// M2[1] = M22 = 2*D22
// M2[2] = M33 = 2*D33
// M2[3] = M12 = 2*D12
// M2[4] = M13 = 2*D13
// M2[5] = M23 = 2*D23

// M4[0]  = M1111 = 12*D11*D11 + 24*D1111
// M4[1]  = M2222 = 12*D22*D22 + 24*D2222
// M4[2]  = M3333 = 12*D33*D33 + 24*D3333

// M4[3]  = M1112 = 12*D11*D12 + 24*D1112
// M4[4]  = M1113 = 12*D11*D13 + 24*D1113
// M4[5]  = M2221 = 12*D22*D12 + 24*D2221
// M4[6]  = M2223 = 12*D22*D23 + 24*D2223
// M4[7]  = M3331 = 12*D33*D13 + 24*D3331
// M4[8]  = M3332 = 12*D33*D23 + 24*D3332

// M4[9]  = M1122 = 4*D11*D22 + 8*D12*D12 + 24*D1122
// M4[10] = M1133 = 4*D11*D33 + 8*D13*D13 + 24*D1133
// M4[11] = M2233 = 4*D22*D33 + 8*D23*D23 + 24*D2233

// M4[12] = M1123 = 4*D11*D23 + 8*D12*D13 + 24*D1123
// M4[13] = M2213 = 4*D22*D13 + 8*D12*D23 + 24*D2213
// M4[14] = M3312 = 4*D33*D12 + 8*D13*D23 + 24*D3312

// see `Technical Notes' for details
// compute all 2nd order (M2) and 4th order (M4) moments 
// using the 2nd and 4th order tensors.
void second_and_fourth_order_moments(float *x, float *y, float *M2, float *M4, int verbose)
{

// these are wrong because of new definitions for x and possibly y.
   M2[0] = 2*x[0];
   M2[1] = 2*x[3];
   M2[2] = 2*x[5];
   M2[3] = 2*x[1];
   M2[4] = 2*x[2];
   M2[5] = 2*x[4];

   M4[0]  = 12*x[0]*x[0] + 24*y[0];
   M4[1]  = 12*x[1]*x[1] + 24*y[1];
   M4[2]  = 12*x[2]*x[2] + 24*y[2];

   M4[3]  = 12*x[0]*x[3] + 24*y[3];
   M4[4]  = 12*x[0]*x[4] + 24*y[4];
   M4[5]  = 12*x[1]*x[3] + 24*y[5];
   M4[6]  = 12*x[1]*x[5] + 24*y[6];
   M4[7]  = 12*x[2]*x[4] + 24*y[7];
   M4[8]  = 12*x[2]*x[5] + 24*y[8];

   M4[9]  = 4*x[0]*x[1] + 8*x[3]*x[3] + 24*y[9];
   M4[10] = 4*x[0]*x[2] + 8*x[4]*x[4] + 24*y[10];
   M4[11] = 4*x[1]*x[2] + 8*x[5]*x[5] + 24*y[11];

   M4[12] = 4*x[0]*x[5] + 8*x[3]*x[4] + 24*y[12];
   M4[13] = 4*x[1]*x[4] + 8*x[3]*x[5] + 24*y[13];
   M4[14] = 4*x[2]*x[3] + 8*x[4]*x[5] + 24*y[14];

if(verbose) printf("%f %f %f\n",x[1], M4[1], y[1]);
/*
*/
   return;
}

void second_and_fourth_order_moments(float *x, float *M2, float *M4)
{
   M2[0] = 2*x[1]; // M11 = 2*D11
   M2[1] = 2*x[4]; // M22 = 2*D22
   M2[2] = 2*x[6]; // M33 = 2*D33
   M2[3] = 2*x[2]; // M12 = 2*D12
   M2[4] = 2*x[3]; // M13 = 2*D13
   M2[5] = 2*x[5]; // M23 = 2*D23

   M4[0]  = 12*x[1]*x[1] + 24*x[7]; // M1111 = 12*D11*D11 + 24*D1111
   M4[1]  = 12*x[4]*x[4] + 24*x[8]; // M2222 = 12*D22*D22 + 24*D2222
   M4[2]  = 12*x[6]*x[6] + 24*x[9]; // M3333 = 12*D33*D33 + 24*D3333

   M4[3]  = 12*x[1]*x[2] + 24*x[10]; // M1112 = 12*D11*D12 + 24*D1112
   M4[4]  = 12*x[1]*x[3] + 24*x[11]; // M1113 = 12*D11*D13 + 24*D1113
   M4[5]  = 12*x[4]*x[2] + 24*x[12]; // M2221 = 12*D22*D12 + 24*D2221
   M4[6]  = 12*x[4]*x[5] + 24*x[13]; // M2223 = 12*D22*D23 + 24*D2223
   M4[7]  = 12*x[6]*x[3] + 24*x[14]; // M3331 = 12*D33*D13 + 24*D3331
   M4[8]  = 12*x[6]*x[5] + 24*x[15]; // M3332 = 12*D33*D23 + 24*D3332

   M4[9]  = 4*x[1]*x[4] + 8*x[2]*x[2] + 24*x[16];  // M1122 = 4*D11*D22 + 8*D12*D12 + 24*D1122
   M4[10] = 4*x[1]*x[6] + 8*x[3]*x[3] + 24*x[17];  // M1133 = 4*D11*D33 + 8*D13*D13 + 24*D1133
   M4[11] = 4*x[4]*x[6] + 8*x[5]*x[5] + 24*x[18];  // M2233 = 4*D22*D33 + 8*D23*D23 + 24*D2233

   M4[12] = 4*x[1]*x[5] + 8*x[2]*x[3] + 24*x[19];  // M1123 = 4*D11*D23 + 8*D12*D13 + 24*D1123
   M4[13] = 4*x[4]*x[3] + 8*x[2]*x[5] + 24*x[20];  // M2213 = 4*D22*D13 + 8*D12*D23 + 24*D2213
   M4[14] = 4*x[6]*x[2] + 8*x[3]*x[5] + 24*x[21];  // M3312 = 4*D33*D12 + 8*D13*D23 + 24*D3312

   return;
}

float mardia_kurtosis(float *x1, float *x2, int v)
{
   float K;
   float M2[6];
   float M4[15];
   float det;
   float S11, S22, S33, S12, S13, S23;

   second_and_fourth_order_moments(x1, x2, M2, M4, v);

   det = M2[0]*M2[1]*M2[2] + 2.0*M2[3]*M2[4]*M2[5] - M2[3]*M2[3]*M2[2] - M2[4]*M2[4]*M2[1] - M2[5]*M2[5]*M2[0];

if(v) printf("det = %f\n",det);

   if(det == 0.0) return(0.0);

   S11 = ( M2[1]*M2[2] - M2[5]*M2[5] )/ det;
   S22 = ( M2[0]*M2[2] - M2[4]*M2[4] )/ det;
   S33 = ( M2[0]*M2[1] - M2[3]*M2[3] )/ det;

   S12 = ( M2[4]*M2[5] - M2[3]*M2[2] )/ det;
   S13 = ( M2[3]*M2[5] - M2[4]*M2[1] )/ det;
   S23 = ( M2[3]*M2[4] - M2[5]*M2[0] )/ det;

if(v) printf("S22=%f M4=%f\n",S22,M4[1]);
if(v) printf("1 %f\n",S22*S22*M4[1]);

   K = S11*S11*M4[0] + S22*S22*M4[1] + S33*S33*M4[2] + 4*S11*S12*M4[3] + 4*S11*S13*M4[4] +
   4*S22*S12*M4[5] + 4*S22*S23*M4[6] + 4*S33*S13*M4[7] + 4*S33*S23*M4[8] +
   (2*S11*S22 + 4*S12*S12)*M4[9] + 
   (2*S11*S33 + 4*S13*S13)*M4[10] +
   (2*S22*S33 + 4*S23*S23)*M4[11] +
   (4*S11*S23 + 8*S12*S13)*M4[12] +
   (4*S22*S13 + 8*S12*S23)*M4[13] +
   (4*S33*S12 + 8*S13*S23)*M4[14] - 15.0;

   return(K);
}

float mardia_kurtosis(float *x)
{
   float K;
   float M2[6];
   float M4[15];
   float det;
   float S11, S22, S33, S12, S13, S23;

   second_and_fourth_order_moments(x, M2, M4);

   det = M2[0]*M2[1]*M2[2] + 2.0*M2[3]*M2[4]*M2[5] - M2[3]*M2[3]*M2[2] - M2[4]*M2[4]*M2[1] - M2[5]*M2[5]*M2[0];

   if(det == 0.0) return(0.0);

   S11 = ( M2[1]*M2[2] - M2[5]*M2[5] )/ det;
   S22 = ( M2[0]*M2[2] - M2[4]*M2[4] )/ det;
   S33 = ( M2[0]*M2[1] - M2[3]*M2[3] )/ det;

   S12 = ( M2[4]*M2[5] - M2[3]*M2[2] )/ det;
   S13 = ( M2[3]*M2[5] - M2[4]*M2[1] )/ det;
   S23 = ( M2[3]*M2[4] - M2[5]*M2[0] )/ det;

   K = S11*S11*M4[0] + S22*S22*M4[1] + S33*S33*M4[2] + 4*S11*S12*M4[3] + 4*S11*S13*M4[4] +
   4*S22*S12*M4[5] + 4*S22*S23*M4[6] + 4*S33*S13*M4[7] + 4*S33*S23*M4[8] +
   (2*S11*S22 + 4*S12*S12)*M4[9] + 
   (2*S11*S33 + 4*S13*S13)*M4[10] +
   (2*S22*S33 + 4*S23*S23)*M4[11] +
   (4*S11*S23 + 8*S12*S13)*M4[12] +
   (4*S22*S13 + 8*S12*S23)*M4[13] +
   (4*S33*S12 + 8*S13*S23)*M4[14] - 15.0;

   return(K);
}
