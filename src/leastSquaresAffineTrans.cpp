#include "babak_lib.h"
#include "bbk_linear_algebra.h"

void leastSquaresAffineTrans(float4 *P, float4 *Q, int4 n, float4 *A)
{
   float4 pavg[3]; // average of columns of P
   float4 qavg[3]; // average of columns of Q
   float4 PPT[9];
   float4 QPT[9];
   float4 *iPPT;
   float4 a[9]; // upper left 3x3 submatrix of A
   float4 d[3];

   // compute pavg and qavg
   avgCol(P, 3, n, pavg);
   avgCol(Q, 3, n, qavg);

   // subtract pavg and qavg from columns of P and Q
   // Note: P and Q are modified
   subtractAvgCol(P, 3, n, pavg);
   subtractAvgCol(Q, 3, n, qavg);

   mat_mat_trans(P,3,n,P,3,PPT);
   mat_mat_trans(Q,3,n,P,3,QPT);
   iPPT = inv3x3(PPT);

   // computes QP'*(PP')^-1
   multi(QPT,3,3,iPPT,3,3,a);

   // computes: d = qavg - a*pavg;
   multi(a,3,3,pavg,3,1,d);
   for(int4 i=0; i<3; i++) d[i] = qavg[i]-d[i];
   
   A[0]=a[0]; A[1]=a[1]; A[2]=a[2]; A[3]=d[0];
   A[4]=a[3]; A[5]=a[4]; A[6]=a[5]; A[7]=d[1];
   A[8]=a[6]; A[9]=a[7]; A[10]=a[8]; A[11]=d[2];
   A[12]=0.0; A[13]=0.0; A[14]=0.0; A[15]=1.0;

   free(iPPT);
}
