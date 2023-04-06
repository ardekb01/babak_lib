#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <babak_lib.h>
#include <assert.h>

#define _matrixCom

void s3adjoint(double *A, double *ADJ);
void s3adjoint(float *A, float *ADJ);

void subtractRowAvg(float *X, int N, int P, float *X0);

void copyVector(float *v1, float *v2, int n);
void subtractVector(float *v1, float *v2, int n);
void normalizeVector(float *x, int n, double *norm);

int centerMatrixRow(float *X, int N, int P, float *avg);
int centerMatrixRow(float *X, int N, int P);

int avgRow(float *X, int N, int P, float *avg);
int avgRow(double *X, int N, int P, double *avg);
int varRow(float *X, int N, int P, float *avg, float *var);
int varRow(double *X, int N, int P, double *avg, double *var);

void transpose_matrix(float *A, int N,  int M);
void transpose_matrix(float *A, int N,  int M, float *AT);

/*****************************************************************************/
/* computes the determinant of the 3x3 matrix A */
double det3(double *A);

/* computes the transpose of the NxM matrix A */
float *trans(float *A, int N,  int M);
double *trans(double *A, int N,  int M);
/*****************************************************************************/

/* computes the determinant of the 4x4 matrix A */
float det3(float *A);
float det4(float *A);
double det4(double *A);

/* computes the inverse of the 4x4 matrix A */
float *inv4(float *A);
double *inv4(double *A);
double *inv3(double *A);
float *inv3(float *A);
void inv3(float *A, float *invA);
float *inv2(float *A);
double *inv2(double *A);

//////////////////////////////////////////////////////////////////

void subtractRowAvg(float *X, int N, int P, float *X0)
{
	for(int i=0; i<N; i++) 
	{
		for(int j=0; j<P; j++) 
			X[i*P +j] -= X0[j];
	}
}

//////////////////////////////////////////////////////////////////

void copyVector(float *v1, float *v2, int n)
{
	for(int i=0; i<n; i++) v1[i] = v2[i];
}

void subtractVector(float *v1, float *v2, int n)
{
	for(int i=0; i<n; i++)
		v1[i] -= v2[i];
}

//////////////////////////////////////////////////////////////////////////////

int subtractAvgCol(float *X, int N, int P, float *avg)
{
   float *rowptr;

   if(N<=0 || P<=0 || X==NULL || avg==NULL) return(1);

   for(int i=0; i<N; i++) 
   {
      rowptr = X+i*P;

      for(int j=0; j<P; j++) 
      {
         rowptr[j] -= avg[i];
      }
   }

   return(0);
}

int centerMatrixRow(float *X, int N, int P, float *avg)
{
	if(N<=0 || P<=0 || X==NULL || avg==NULL) return(1);

	for(int i=0; i<N; i++) 
	{
		for(int j=0; j<P; j++) 
			X[i*P +j] -= avg[j];
	}

	return(0);
}

int centerMatrixRow(float *X, int N, int P)
{
	float *avg;

	if(N<=0 || P<=0 || X==NULL) return(1);

	avg = (float *)calloc(P, sizeof(float));

	avgRow(X, N, P, avg);

	for(int i=0; i<N; i++) 
	{
		for(int j=0; j<P; j++) 
			X[i*P +j] -= avg[j];
	}

	free(avg);

	return(0);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Returns 1 if an error condition occurs, 0 otherwise
int avgRow(double *X, int N, int P, double *avg, char *rowmask)
{
	int n=0;

	// j is the column index
	// i is the row index
	// P is the number of columns
	// N is the number of rows 

	if(N<=0 || P<=0 || X==NULL || avg==NULL || rowmask==NULL) return(1);

	for(int i=0; i<N; i++) if(rowmask[i]==0) n++;

	if(n==0) return(1);

	for(int j=0; j<P; j++) 		avg[j] = 0.0;

	for(int i=0; i<N; i++) 		// add all rows
	if( rowmask[i] == 0 )
	{
		for(int j=0; j<P; j++) 		// for each column
			avg[j] += X[i*P + j];
	}

	for(int j=0; j<P; j++) 	avg[j] /= n;

	return(0);
}

// Returns 1 if an error condition occurs, 0 otherwise
int avgRow(float *X, int N, int P, float *avg, char *rowmask)
{
	int n=0;

	// j is the column index
	// i is the row index
	// P is the number of columns
	// N is the number of rows 

	if(N<=0 || P<=0 || X==NULL || avg==NULL || rowmask==NULL) return(1);

	for(int i=0; i<N; i++) if(rowmask[i]==0) n++;

	if(n==0) return(1);

	for(int j=0; j<P; j++) 		avg[j] = 0.0;

	for(int i=0; i<N; i++) 		// add all rows
	if( rowmask[i] == 0 )
	{
		for(int j=0; j<P; j++) 		// for each column
			avg[j] += X[i*P + j];
	}

	for(int j=0; j<P; j++) 	avg[j] /= n;

	return(0);
}

// Returns 1 if an error condition occurs, 0 otherwise
int avgRow(float *X, int N, int P, float *avg)
{
   // j is the column index
   // i is the row index
   // P is the number of columns
   // N is the number of rows 

   if(N<=0 || P<=0 || X==NULL || avg==NULL) return(1);

   for(int j=0; j<P; j++) 		// for each column
   {
      avg[j]=0.0;

      for(int i=0; i<N; i++) 		// add all rows
      {
         avg[j] += X[i*P + j];
      }

      avg[j] /= N;
   }

   return(0);
}

// returns an N-vector in "avg" containing the average of the P columns of X
// X: N*P matrix
// avg: N vector
// Averages X across the P dimension
int avgCol(float *X, int N, int P, float *avg)
{
   float *rowptr; // points to the begining of a row
   float rowsum;

   // j is the column index
   // i is the row index
   // P is the number of columns
   // N is the number of rows 

   if(N<=0 || P<=0 || X==NULL || avg==NULL) return(1);

   for(int i=0; i<N; i++) // for each row
   {
      rowsum=0.0;
      rowptr = X+i*P;

      for(int j=0; j<P; j++) // add all columns
      {
         rowsum += rowptr[j];
      }

      avg[i] = rowsum/P;
   }

   return(0);
}

// Returns 1 if an error condition occurs, 0 otherwise
int avgRow(double *X, int N, int P, double *avg)
{
	// j is the column index
	// i is the row index
	// P is the number of columns
	// N is the number of rows 

	if(N<=0 || P<=0 || X==NULL || avg==NULL) return(1);

	for(int j=0; j<P; j++) 		// for each column
	{
		avg[j]=0.0;
		for(int i=0; i<N; i++) 		// add all rows
			avg[j] += X[i*P + j];
		avg[j] /= N;
	}

	return(0);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Returns 1 if an error condition occurs, 0 otherwise
int varRow(double *X, int N, int P, double *avg, double *var, char *rowmask)
{
	int n=0;

	// j is the column index
	// i is the row index
	// P is the number of columns
	// N is the number of rows 

	if(N<2 || P<1 || X==NULL || avg==NULL || var==NULL || rowmask==NULL) return(1);

	for(int i=0; i<N; i++) if(rowmask[i]==0) n++;

	if( n<2 ) return(1);

	for(int j=0; j<P; j++) 	var[j] = 0.0;

	for(int i=0; i<N; i++) 		// add all rows
	if(rowmask[i]==0)
	{
		for(int j=0; j<P; j++) 		// for each column
			var[j] += (X[i*P + j]-avg[j])*(X[i*P + j]-avg[j]);
	}

	for(int j=0; j<P; j++) var[j] /= (n-1);

	return(0);
}

// Returns 1 if an error condition occurs, 0 otherwise
int varRow(float *X, int N, int P, float *avg, float *var, char *rowmask)
{
	int n=0;

	// j is the column index
	// i is the row index
	// P is the number of columns
	// N is the number of rows 

	if(N<2 || P<1 || X==NULL || avg==NULL || var==NULL || rowmask==NULL) return(1);

	for(int i=0; i<N; i++) if(rowmask[i]==0) n++;

	if( n<2 ) return(1);

	for(int j=0; j<P; j++) 	var[j] = 0.0;

	for(int i=0; i<N; i++) 		// add all rows
	if(rowmask[i]==0)
	{
		for(int j=0; j<P; j++) 		// for each column
			var[j] += (X[i*P + j]-avg[j])*(X[i*P + j]-avg[j]);
	}

	for(int j=0; j<P; j++) var[j] /= (n-1);

	return(0);
}

// Returns 1 if an error condition occurs, 0 otherwise
int ssdRow(float *X, int N, int P, float *avg, float *ssd, char *rowmask)
{
	int n=0;

	// j is the column index
	// i is the row index
	// P is the number of columns
	// N is the number of rows 

	if(N<2 || P<1 || X==NULL || avg==NULL || ssd==NULL || rowmask==NULL) return(1);

	for(int i=0; i<N; i++) if(rowmask[i]==0) n++;

	if( n<2 ) return(1);

	for(int j=0; j<P; j++) 	ssd[j] = 0.0;

	for(int i=0; i<N; i++) 		// add all rows
	if(rowmask[i]==0)
	{
		for(int j=0; j<P; j++) 		// for each column
			ssd[j] += (X[i*P + j]-avg[j])*(X[i*P + j]-avg[j]);
	}

	return(0);
}

// Returns 1 if an error condition occurs, 0 otherwise
int varRow(float *X, int N, int P, float *avg, float *var)
{
	// j is the column index
	// i is the row index
	// P is the number of columns
	// N is the number of rows 

	if(N<=1 || P<=0 || X==NULL || avg==NULL || var==NULL) return(1);

	for(int j=0; j<P; j++) 		// for each column
	{
		var[j]=0.0;

		for(int i=0; i<N; i++) 		// add all rows
			var[j] += (X[i*P + j]-avg[j])*(X[i*P + j]-avg[j]);

		var[j] /= (N-1);
	}

	return(0);
}

// Returns 1 if an error condition occurs, 0 otherwise
int ssdRow(float *X, int N, int P, float *avg, float *ssd)
{
	if(N<=1 || P<=0 || X==NULL || avg==NULL || ssd==NULL) return(1);

	for(int j=0; j<P; j++) 		// for each column
	{
		ssd[j]=0.0;

		for(int i=0; i<N; i++) 		// add all rows
			ssd[j] += (X[i*P + j]-avg[j])*(X[i*P + j]-avg[j]);
	}

	return(0);
}

// Returns 1 if an error condition occurs, 0 otherwise
int varRow(double *X, int N, int P, double *avg, double *var)
{
	// j is the column index
	// i is the row index
	// P is the number of columns
	// N is the number of rows 

	if(N<=1 || P<=0 || X==NULL || avg==NULL || var==NULL) return(1);

	for(int j=0; j<P; j++) 		// for each column
	{
		var[j]=0.0;

		for(int i=0; i<N; i++) 		// add all rows
			var[j] += (X[i*P + j]-avg[j])*(X[i*P + j]-avg[j]);

		var[j] /= (N-1);
	}

	return(0);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// A is a pxp matrix
// x is a px1 matrix
// The function retuns xt * A * x where xt denotes the transpose of x.
double xtAx(float *A, double *x, int p)
{
	double sum=0;

	for(int i=0; i<p; i++)
	for(int j=0; j<p; j++)
		sum += A[i*p + j]*x[i]*x[j];

	return(sum);
}
//////////////////////////////////////////////////////////////////////////////


// Multiplies floating type matrices A and B. The result is returned in C=AB.
//The number of columns in A must be equal to the number of rows in B
//(i.e., jA=iB), otherwise error results.  
void multi(float *A,int iA,int jA,float *B,int iB,int jB,float *C);
void multi(double *A,int iA,int jA,double *B,int iB,int jB,double *C);
/*****************************************************************************/

double vectorNorm(float *x, int n)
{
   double norm=0.0;

   if(x==NULL)
   {
      printf("\n\nWarning: NULL array passed to vectorNorm() function.\n\n");
      return(0.0);
   }

   if( n<=0 )
   {
      printf("\n\nWarning: Non-positive array dimension argument passed to vectorNorm() function.\n\n");
      return(0.0);
   }

   for(int i=0; i<n; i++)
      norm += x[i]*x[i];

   norm = sqrt(norm);

   return(norm);
}

void normalizeVector(float *x, int n)
{
   double norm=0.0;

   if(x==NULL)
   {
      printf("\n\nWarning: NULL array passed to normalizeVector() function.\n\n");
      return;
   }

   if( n<=0 )
   {
      printf("\n\nWarning: Non-positive array dimension argument passed to normalizeVector() function.\n\n");
      return;
   }

   norm=vectorNorm(x, n);

   if( norm == 0.0) return;

   for(int i=0; i<n; i++)
   {
      x[i] = (float)(x[i]/norm);
   }

   return;
}

void normalizeVector(float *x, int n, double *norm)
{
		*norm=vectorNorm(x, n);

        if( *norm == 0.0) return;

        for(int i=0; i<n; i++)
                x[i] = (float)( x[i] / (*norm) );

        return;
}

///////////////////////////////////////////////////////////////
float normalize(float *s, int n)
{
        float norm=0.0;

        for(int i=0; i<n; i++)
                norm += s[i]*s[i];

        norm = sqrtf(norm);

        if( norm < ESMALL) return(0.0);

        for(int i=0; i<n; i++)
                s[i] = s[i]/norm;

        return(norm);
}

double normalize(double *s, int n)
{
   double norm=0.0;

   for(int i=0; i<n; i++)
      norm += s[i]*s[i];

   norm = sqrt(norm);

   if( norm < ESMALL) return(0.0);

   for(int i=0; i<n; i++)
      s[i] /= norm;

   return(norm);
}
///////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// M: Input 3x3 symmetric matrix with 6 unique elements.
// V: Output 6x1 vector containing the unique elements of M.
// Memory for V should be allocated before calling this function.
// Implements one of several ways that V can be defined.
// V should be considered as a vector form of the symmetric matrix M.
////////////////////////////////////////////////////////////////////////////////
void s3mat_to_vec(float *M, float *V)
{
   V[0] = M[0];
   V[1] = M[4];
   V[2] = M[8];

   V[3] = M[3];
   V[4] = M[7];

   V[5] = M[6];
   
   return;
}

void s3mat_to_vec(float *M, double *V)
{
   V[0] = M[0];
   V[1] = M[4];
   V[2] = M[8];

   V[3] = M[3];
   V[4] = M[7];

   V[5] = M[6];
   
   return;
}

void s3mat_to_vec(double *M, double *V)
{
   V[0] = M[0];
   V[1] = M[4];
   V[2] = M[8];

   V[3] = M[3];
   V[4] = M[7];

   V[5] = M[6];
   
   return;
}

void s3vec_to_mat(double *M, double *V)
{
   M[0] = V[0];
   M[4] = V[1];
   M[8] = V[2];
   M[3] = V[3];
   M[7] = V[4];
   M[6] = V[5];
   M[1] = V[3];
   M[2] = V[5];
   M[5] = V[4];
   
   return;
}

void s3vec_to_mat(float *M, float *V)
{
   M[0] = V[0];
   M[4] = V[1];
   M[8] = V[2];
   M[3] = V[3];
   M[7] = V[4];
   M[6] = V[5];
   M[1] = V[3];
   M[2] = V[5];
   M[5] = V[4];
   
   return;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Computes |A| where A is a 3x3 symmetric matrix in vector form
//      A0 A3 A5
// A =  A3 A1 A4
//      A5 A4 A2
// Given as: A = [A0 A1 A2 A3 A4 A5]
////////////////////////////////////////////////////////////////////////////////////////////////
double s3det(double *A)
{
   return (A[0]*A[1]*A[2] + 2.0*A[3]*A[4]*A[5] 
   - A[0]*A[4]*A[4] 
   - A[3]*A[3]*A[2] 
   - A[5]*A[1]*A[5] );
}

float s3det(float *A)
{
   return( A[0]*A[1]*A[2] + 2.0*A[3]*A[4]*A[5] 
   - A[0]*A[4]*A[4] 
   - A[3]*A[3]*A[2] 
   - A[5]*A[1]*A[5] );
}
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Computes d|A|/dA where A is a 3x3 symmetric matrix in vector form
//      A0 A3 A5
// A =  A3 A1 A4
//      A5 A4 A2
// Given as: A = [A0 A1 A2 A3 A4 A5]
// Output B=d|A|/dA given in 6x1 vector form
////////////////////////////////////////////////////////////////////////////////////////////////
void ds3det(double *A, double *B)
{
   B[0] = A[1]*A[2] - A[4]*A[4];
   B[1] = A[0]*A[2] - A[5]*A[5];
   B[2] = A[0]*A[1] - A[3]*A[3];
   B[3] = 2.0*(A[4]*A[5] - A[2]*A[3]);
   B[4] = 2.0*(A[3]*A[5] - A[0]*A[4]);
   B[5] = 2.0*(A[3]*A[4] - A[1]*A[5]);

   return;
}

void ds3det(float *A, float *B)
{
   B[0] = A[1]*A[2] - A[4]*A[4];
   B[1] = A[0]*A[2] - A[5]*A[5];
   B[2] = A[0]*A[1] - A[3]*A[3];
   B[3] = 2.0*(A[4]*A[5] - A[2]*A[3]);
   B[4] = 2.0*(A[3]*A[5] - A[0]*A[4]);
   B[5] = 2.0*(A[3]*A[4] - A[1]*A[5]);

   return;
}
////////////////////////////////////////////////////////////////////////////////////////////////

// Computes tr[A*B]
// The trace of A*B
// A and B are in S(3) given as 6x1 vectors
double s3tr(double *A, double *B)
{
   return(
      A[0]*B[0] +
      A[1]*B[1] +
      A[2]*B[2] +
      2.0*A[3]*B[3] +
      2.0*A[4]*B[4] +
      2.0*A[5]*B[5]
   );
}

float s3tr(float *A, float *B)
{
   return(
      A[0]*B[0] +
      A[1]*B[1] +
      A[2]*B[2] +
      2.0*A[3]*B[3] +
      2.0*A[4]*B[4] +
      2.0*A[5]*B[5]
   );
}

// Given matrices A and B as inputs, this function computes dtr( adj(A)*B)/dA
// A and B are S(3) given as 6x1 vectors
void compute_dc1(double *A, double *B, double *C)
{
   C[0] = B[1]*A[2] + B[2]*A[1] - 2.0*B[4]*A[4];
   C[1] = B[0]*A[2] + B[2]*A[0] - 2.0*B[5]*A[5];
   C[2] = B[0]*A[1] + B[1]*A[0] - 2.0*B[3]*A[3];
   C[3] = 2.0*( -B[2]*A[3] - B[3]*A[2] + B[4]*A[5] + B[5]*A[4] );
   C[4] = 2.0*( -B[0]*A[4] + B[3]*A[5] - B[4]*A[0] + B[5]*A[3] );
   C[5] = 2.0*( -B[1]*A[5] + B[3]*A[4] + B[4]*A[3] - B[5]*A[1] );

   return;
}

////////////////////////////////////////////////////////////////////////////////////////////////

// A and B are 6x1 vectors in S(3)
// B=adjoint(A)
void s3adjoint(double *A, double *B)
{
   B[0] = A[1]*A[2] - A[4]*A[4];
   B[1] = A[0]*A[2] - A[5]*A[5];
   B[2] = A[0]*A[1] - A[3]*A[3];
   B[3] = A[4]*A[5] - A[2]*A[3];
   B[4] = A[3]*A[5] - A[0]*A[4];
   B[5] = A[3]*A[4] - A[1]*A[5];

   return;
}

void s3adjoint(float *A, float *B)
{
   B[0] = A[1]*A[2] - A[4]*A[4];
   B[1] = A[0]*A[2] - A[5]*A[5];
   B[2] = A[0]*A[1] - A[3]*A[3];
   B[3] = A[4]*A[5] - A[2]*A[3];
   B[4] = A[3]*A[5] - A[0]*A[4];
   B[5] = A[3]*A[4] - A[1]*A[5];

   return;
}

//////////////////////////////////////////////////////
void s3inv(double *A, double *invA)
{
   double det;

   det = s3det(A);

   s3adjoint(A, invA);

   if(det != 0.0)
   {
      for(int i=0; i<6; i++)
      {
         invA[i]/=det;
      }
   }
   else
   {
      for(int i=0; i<6; i++)
      {
         invA[i]=0.0;
      }
   }
}

int s3inv(float *A, float *invA)
{
   double det;

   det = s3det(A);

   s3adjoint(A, invA);

   if(det != 0.0)
   {
      for(int i=0; i<6; i++)
      {
         invA[i]/=det;
      }
      return 1;
   }
   else
   {
      for(int i=0; i<6; i++)
      {
         invA[i]=0.0;
      }
      return 0;
   }
}
//////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
// D and F 3x3 symmetric matrix as 6x1 vector
// L 3x1 vector of eigenvalues
// dL 6x1 matrix of derivatives wrt D
// c3*L^3 - c2*L^2 + c1*L - c0 = 0
// c3=|F|
// c2=trace( adjoint(F)*D )
// c1=trace( adjoint(D)*F )
// c0=|D|
///////////////////////////////////////////////////////////////
void differentiate_distance(double *D, double *F, double *L, double *dL)
{
   double c3=0.0;
   double c2=0.0;
   double c1=0.0;

   double dc2[6]; // derivative of c2 wrt matrix D
   double dc1[6]; // derivative of c1 wrt matrix D
   double dc0[6]; // derivative of c0 wrt matrix D

   double adjF[6]; // adjoint of matrix F
   double adjD[6]; // adjoint of matrix D

   double dum, l, l2, lnl;

   s3adjoint(F, adjF);
   s3adjoint(D, adjD);

   c3 = s3det(F);
   c2 = s3tr(adjF, D);
   c1 = s3tr(adjD, F);

   // compute dc0/dD;
   for(int i=0; i<3; i++) dc0[i]=adjD[i];
   for(int i=3; i<6; i++) dc0[i]=2.0*adjD[i];

   // compute dc1/dD
   compute_dc1(D, F, dc1);

   // compute dc2/dD
   for(int i=0; i<3; i++) dc2[i]=adjF[i];
   for(int i=3; i<6; i++) dc2[i]=2.0*adjF[i];

   for(int i=0; i<6; i++)
   {
      dL[i] = 0.0;
   }

   for(int n=0; n<3; n++)
   {
      
      l = L[n];
      l2 = l*l;
      lnl = log(l);

      dum =    (3.0*c3*l2 - 2.0*c2*l + c1)*l/lnl;
      // dum = (3.0*c3*l2 - 2.0*c2*l + c1);

      if(dum!=0.0)
      {
         for(int i=0; i<6; i++)
         {
            dL[i]   += (l2*dc2[i] - l*dc1[i] + dc0[i])/dum;
            // dL[i] = (l2*dc2[i] - l*dc1[i] + dc0[i])/dum;
         }
      }

      //printf("\n");
      //for(int i=0; i<6; i++)
      //   printf("%lf\n",dL[i]);
   }
   //printf("\n");
   //for(int i=0; i<6; i++)
   //   printf("%lf\n",dL[i]);
}

///////////////////////////////////////////////////////////////
// compute the eigenvalues of the symmetric matrix A
void s3eigenval(double *A, double *L)
{
   double p;
   double q;
   double B[6];
   double r;
   double phi;

   static double PI=4.0*atan(1.0);

   p = A[3]*A[3] + A[5]*A[5] + A[4]*A[4];

   if( p == 0.0 )
   {
      // Covers all six ranking possibilities so that eigenvalues 
      // are reported in increasing order.
      if(A[0]>A[1] && A[1]>A[2]) // A0>A1>A2
      {
         L[0] = A[0];   
         L[1] = A[1];   
         L[2] = A[2];   
      }
      else if(A[0]>A[2] && A[2]>A[1]) // A0>A2>A1
      {
         L[0] = A[0];   
         L[1] = A[2];   
         L[2] = A[1];   
      }
      else if(A[1]>A[0] && A[0]>A[2]) // A1>A0>A2
      {
         L[0] = A[1];   
         L[1] = A[0];   
         L[2] = A[2];   
      }
      else if(A[1]>A[2] && A[2]>A[0]) // A1>A2>A0
      {
         L[0] = A[1];   
         L[1] = A[2];   
         L[2] = A[0];   
      }
      else if(A[2]>A[0] && A[0]>A[1]) // A2>A0>A1
      {
         L[0] = A[2];   
         L[1] = A[0];   
         L[2] = A[1];   
      }
      else if(A[2]>A[1] && A[1]>A[0]) // A2>A1>A0
      {
         L[0] = A[2];   
         L[1] = A[1];   
         L[2] = A[0];   
      }
   }
   else
   {
      q = (A[0] + A[1] + A[2])/3.0;

      p = (A[0] - q)*(A[0] - q) + (A[1] - q)*(A[1] - q) + (A[2] - q)*(A[2] - q) + 2.0*p;
      p = sqrt(p/6.0);

      B[0] = (A[0] - q)/p;
      B[1] = (A[1] - q)/p;
      B[2] = (A[2] - q)/p;
      B[3] = A[3]/p;
      B[4] = A[4]/p;
      B[5] = A[5]/p;

      r = s3det(B)/2.0;

      if( r <= -1.0)
      {
         phi = PI/3.0;
      }
      else if( r >= 1.0)
      {
         phi = 0.0;
      }
      else
      {
         phi = acos(r)/3.0;
      }

      L[0] = q + 2.0*p*cos(phi);
      L[2] = q + 2.0*p*cos(phi + PI*2.0/3.0);
      L[1] = 3.0*q - L[0] - L[2];
   }

   return;
}

void s3eigenval(float *A, float *L)
{
   float p;
   float q;
   float B[6];
   float r;
   float phi;

   static float PI=4.0*atanf(1.0);

   p = A[3]*A[3] + A[5]*A[5] + A[4]*A[4];

   if( p == 0.0 )
   {
      // Covers all six ranking possibilities so that eigenvalues 
      // are reported in increasing order.
      if(A[0]>A[1] && A[1]>A[2]) // A0>A1>A2
      {
         L[0] = A[0];   
         L[1] = A[1];   
         L[2] = A[2];   
      }
      else if(A[0]>A[2] && A[2]>A[1]) // A0>A2>A1
      {
         L[0] = A[0];   
         L[1] = A[2];   
         L[2] = A[1];   
      }
      else if(A[1]>A[0] && A[0]>A[2]) // A1>A0>A2
      {
         L[0] = A[1];   
         L[1] = A[0];   
         L[2] = A[2];   
      }
      else if(A[1]>A[2] && A[2]>A[0]) // A1>A2>A0
      {
         L[0] = A[1];   
         L[1] = A[2];   
         L[2] = A[0];   
      }
      else if(A[2]>A[0] && A[0]>A[1]) // A2>A0>A1
      {
         L[0] = A[2];   
         L[1] = A[0];   
         L[2] = A[1];   
      }
      else if(A[2]>A[1] && A[1]>A[0]) // A2>A1>A0
      {
         L[0] = A[2];   
         L[1] = A[1];   
         L[2] = A[0];   
      }
   }
   else
   {
      q = (A[0] + A[1] + A[2])/3.0;

      p = (A[0] - q)*(A[0] - q) + (A[1] - q)*(A[1] - q) + (A[2] - q)*(A[2] - q) + 2.0*p;
      p = sqrtf(p/6.0);

      B[0] = (A[0] - q)/p;
      B[1] = (A[1] - q)/p;
      B[2] = (A[2] - q)/p;
      B[3] = A[3]/p;
      B[4] = A[4]/p;
      B[5] = A[5]/p;

      r = s3det(B)/2.0;

      if( r <= -1.0)
      {
         phi = PI/3.0;
      }
      else if( r >= 1.0)
      {
         phi = 0.0;
      }
      else
      {
         phi = acosf(r)/3.0;
      }

      L[0] = q + 2.0*p*cosf(phi);
      L[2] = q + 2.0*p*cosf(phi + PI*2.0/3.0);
      L[1] = 3.0*q - L[0] - L[2];
   }

   return;
}
///////////////////////////////////////////////////////////////

float det3(float *A)
{
	float det;

	det = (A[0]*A[4]*A[8] + A[1]*A[5]*A[6] + A[2]*A[3]*A[7] - A[2]*A[4]*A[6] - A[0]*A[5]*A[7] - A[1]*A[3]*A[8] );

   return(det);
}

double det3(double *A)
{
	double det;

	det = (A[0]*A[4]*A[8] + A[1]*A[5]*A[6] + A[2]*A[3]*A[7] -
   		 A[2]*A[4]*A[6] - A[0]*A[5]*A[7] - A[1]*A[3]*A[8] );

	return(det);
}

float det4(float *A)
{
   float det;
   float B[9];   
   
   det=0.0;

   B[0]=A[5]; B[1]=A[6]; B[2]=A[7]; 
   B[3]=A[9]; B[4]=A[10]; B[5]=A[11]; 
   B[6]=A[13]; B[7]=A[14]; B[8]=A[15]; 
 
   det += A[0]*det3(B);
   
   B[0]=A[4]; B[1]=A[6]; B[2]=A[7]; 
   B[3]=A[8]; B[4]=A[10]; B[5]=A[11]; 
   B[6]=A[12]; B[7]=A[14]; B[8]=A[15]; 
   
   det -= A[1]*det3(B);

   B[0]=A[4]; B[1]=A[5]; B[2]=A[7]; 
   B[3]=A[8]; B[4]=A[9]; B[5]=A[11]; 
   B[6]=A[12]; B[7]=A[13]; B[8]=A[15]; 
   
   det += A[2]*det3(B);

   B[0]=A[4]; B[1]=A[5]; B[2]=A[6]; 
   B[3]=A[8]; B[4]=A[9]; B[5]=A[10]; 
   B[6]=A[12]; B[7]=A[13]; B[8]=A[14]; 
   
   det -= A[3]*det3(B);

   return(det);
}

double det4(double *A)
{
   double det;
   double B[9];   
   
   det=0.0;

   B[0]=A[5]; B[1]=A[6]; B[2]=A[7]; 
   B[3]=A[9]; B[4]=A[10]; B[5]=A[11]; 
   B[6]=A[13]; B[7]=A[14]; B[8]=A[15]; 
 
   det += A[0]*det3(B);
   
   B[0]=A[4]; B[1]=A[6]; B[2]=A[7]; 
   B[3]=A[8]; B[4]=A[10]; B[5]=A[11]; 
   B[6]=A[12]; B[7]=A[14]; B[8]=A[15]; 
   
   det -= A[1]*det3(B);

   B[0]=A[4]; B[1]=A[5]; B[2]=A[7]; 
   B[3]=A[8]; B[4]=A[9]; B[5]=A[11]; 
   B[6]=A[12]; B[7]=A[13]; B[8]=A[15]; 
   
   det += A[2]*det3(B);

   B[0]=A[4]; B[1]=A[5]; B[2]=A[6]; 
   B[3]=A[8]; B[4]=A[9]; B[5]=A[10]; 
   B[6]=A[12]; B[7]=A[13]; B[8]=A[14]; 
   
   det -= A[3]*det3(B);

   return(det);
}

void transpose_matrix(float *A, int N,  int M)
{
   float *transA;

   transA=(float *)calloc(N*M,sizeof(float));

   for(int i=0;i<N;i++)
   for(int j=0;j<M;j++)
      transA[j*N+i]=A[i*M+j];

   for(int i=0;i<N*M;i++) A[i]=transA[i];

   delete transA;
}

void transpose_matrix(float *A, int N,  int M, float *AT)
{
   for(int i=0;i<N;i++)
   {
      for(int j=0;j<M;j++)
      {
         AT[j*N+i]=A[i*M+j];
      }
   }

   return;
}

float *trans(float *A, int N,  int M)
{
   int i,j;

   float *transA;

   transA=(float *)calloc(N*M,sizeof(float));

   for(i=0;i<N;i++)
   for(j=0;j<M;j++)
      transA[j*N+i]=A[i*M+j];

   return(transA);
}

double *trans(double *A, int N,  int M)
{
   int i,j;

   double *transA;

   transA=(double *)calloc(N*M,sizeof(double));

   for(i=0;i<N;i++)
   for(j=0;j<M;j++)
      transA[j*N+i]=A[i*M+j];

   return(transA);
}

// Note: A is preserved
float *inv4(float *A)
{
   int i;

   float detA;
   float *invA;
   float C[16];
   float B[9];
   float *transC;

   invA=(float *)calloc(16,sizeof(float));

   B[0]=A[5];  B[1]=A[6];  B[2]=A[7]; 
   B[3]=A[9];  B[4]=A[10]; B[5]=A[11]; 
   B[6]=A[13]; B[7]=A[14]; B[8]=A[15]; 
   C[0]=det3(B);

   B[0]=A[4];  B[1]=A[6];  B[2]=A[7]; 
   B[3]=A[8];  B[4]=A[10]; B[5]=A[11]; 
   B[6]=A[12]; B[7]=A[14]; B[8]=A[15]; 
   C[1]=-det3(B);

   B[0]=A[4];  B[1]=A[5];   B[2]=A[7]; 
   B[3]=A[8];  B[4]=A[9];  B[5]=A[11]; 
   B[6]=A[12]; B[7]=A[13]; B[8]=A[15]; 
   C[2]=det3(B);

   B[0]=A[4];  B[1]=A[5];   B[2]=A[6]; 
   B[3]=A[8];  B[4]=A[9];  B[5]=A[10]; 
   B[6]=A[12]; B[7]=A[13]; B[8]=A[14]; 
   C[3]=-det3(B);

   B[0]=A[1];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[9];  B[4]=A[10]; B[5]=A[11];
   B[6]=A[13]; B[7]=A[14]; B[8]=A[15];
   C[4]=-det3(B);

   B[0]=A[0];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[8];  B[4]=A[10]; B[5]=A[11];
   B[6]=A[12]; B[7]=A[14]; B[8]=A[15];
   C[5]=det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[3];
   B[3]=A[8];  B[4]=A[9]; B[5]=A[11];
   B[6]=A[12]; B[7]=A[13]; B[8]=A[15];
   C[6]=-det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[2];
   B[3]=A[8];  B[4]=A[9];  B[5]=A[10];
   B[6]=A[12]; B[7]=A[13]; B[8]=A[14];
   C[7]=det3(B);

   B[0]=A[1];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[5];  B[4]=A[6];  B[5]=A[7];
   B[6]=A[13];  B[7]=A[14];  B[8]=A[15];
   C[8]=det3(B);

   B[0]=A[0];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[4];  B[4]=A[6];  B[5]=A[7];
   B[6]=A[12];  B[7]=A[14];  B[8]=A[15];
   C[9]=-det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[3];
   B[3]=A[4];  B[4]=A[5];  B[5]=A[7];
   B[6]=A[12];  B[7]=A[13];  B[8]=A[15];
   C[10]=det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[2];
   B[3]=A[4];  B[4]=A[5];  B[5]=A[6];
   B[6]=A[12];  B[7]=A[13];  B[8]=A[14];
   C[11]=-det3(B);

   B[0]=A[1];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[5];  B[4]=A[6];  B[5]=A[7];
   B[6]=A[9]; B[7]=A[10]; B[8]=A[11];
   C[12]=-det3(B);

   B[0]=A[0];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[4];  B[4]=A[6];  B[5]=A[7];
   B[6]=A[8]; B[7]=A[10]; B[8]=A[11];
   C[13]=det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[3];
   B[3]=A[4];  B[4]=A[5];  B[5]=A[7];
   B[6]=A[8]; B[7]=A[9]; B[8]=A[11];
   C[14]=-det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[2];
   B[3]=A[4];  B[4]=A[5];  B[5]=A[6];
   B[6]=A[8]; B[7]=A[9]; B[8]=A[10];
   C[15]=det3(B);

   transC=trans(C,4,4);

   detA=det4(A);

   if(detA!=0.0)
   {
      for(i=0;i<16;i++)  
        invA[i]=transC[i]/detA;
   }

   free(transC);

   return(invA);
}

void multi(double *A, int iA, int jA,  float *B, int iB,  int jB,  float *C)
{
	int i,j,k;
	float sum;
	float *D;

	/* check to see if the number of columns in A is equal to the number of rows
	in B, if not return an error message and exit */
	if(jA != iB) {
		printf("Warning: incompatible matrix dimensions\n");
		return;
	}

	D=(float *)calloc(iA*jB,sizeof(float));

	for(i=0;i<iA;i++)
	for(j=0;j<jB;j++) {
		sum=0.0;
		for(k=0;k<iB;k++)      /* or jA */
			sum += A[i*jA+k]*B[k*jB+j];
		D[i*jB+j]=sum;
	}

	for(i=0;i<iA;i++)
	for(j=0;j<jB;j++)
		C[i*jB+j]=D[i*jB+j];

	free(D);
}

void multi(float *A, int iA, int jA,  float *B, int iB,  int jB,  float *C)
{
	int i,j,k;
	float sum;
	float *D;

	/* check to see if the number of columns in A is equal to the number of rows
	in B, if not return an error message and exit */
	if(jA != iB) {
		printf("Warning: incompatible matrix dimensions\n");
		return;
	}

	D=(float *)calloc(iA*jB,sizeof(float));

	for(i=0;i<iA;i++)
	for(j=0;j<jB;j++) {
		sum=0.0;
		for(k=0;k<iB;k++)      /* or jA */
			sum += A[i*jA+k]*B[k*jB+j];
		D[i*jB+j]=sum;
	}

	for(i=0;i<iA;i++)
	for(j=0;j<jB;j++)
		C[i*jB+j]=D[i*jB+j];

	free(D);
}

void multi(double *A, int iA, int jA,  double *B, int iB,  int jB,  double *C)
{
	int i,j,k;
	double sum;
	double *D;

	/* check to see if the number of columns in A is equal to the number of rows
	in B, if not return an error message and exit */
	if(jA != iB) {
		printf("Warning: incompatible matrix dimensions\n");
		return;
	}

	D=(double *)calloc(iA*jB,sizeof(double));

	for(i=0;i<iA;i++)
	for(j=0;j<jB;j++) {
		sum=0.0;
		for(k=0;k<iB;k++)      /* or jA */
			sum += A[i*jA+k]*B[k*jB+j];
	D[i*jB+j]=sum;
	}

	for(i=0;i<iA;i++)
	for(j=0;j<jB;j++)
		C[i*jB+j]=D[i*jB+j];

	free(D);
}

void multi(float *A, int iA, int jA,  double *B, int iB,  int jB,  double *C)
{
	int i,j,k;
	double sum;
	double *D;

	/* check to see if the number of columns in A is equal to the number of rows
	in B, if not return an error message and exit */
	if(jA != iB) {
		printf("Warning: incompatible matrix dimensions\n");
		return;
	}

	D=(double *)calloc(iA*jB,sizeof(double));

	for(i=0;i<iA;i++)
	for(j=0;j<jB;j++) {
		sum=0.0;
		for(k=0;k<iB;k++)      /* or jA */
			sum += A[i*jA+k]*B[k*jB+j];
	D[i*jB+j]=sum;
	}

	for(i=0;i<iA;i++)
	for(j=0;j<jB;j++)
		C[i*jB+j]=D[i*jB+j];

	free(D);
}

double *inv3(double *A)
{
   double *invA;
   double CT[9];
   double detA;

   invA=(double *)calloc(9,sizeof(double));

   CT[0] = A[4]*A[8] - A[5]*A[7];
   CT[1] = A[2]*A[7] - A[1]*A[8];
   CT[2] = A[1]*A[5] - A[2]*A[4];
   CT[3] = A[5]*A[6] - A[3]*A[8];
   CT[4] = A[0]*A[8] - A[2]*A[6];
   CT[5] = A[2]*A[3] - A[0]*A[5];
   CT[6] = A[3]*A[7] - A[4]*A[6];
   CT[7] = A[1]*A[6] - A[0]*A[7];
   CT[8] = A[0]*A[4] - A[1]*A[3];

   detA=det3(A);

   if(detA!=0.0)
   {
      for(int i=0; i<9; i++)  
         invA[i]=CT[i]/detA;
   }

   return(invA);
}

float *inv3(float *A)
{
   float *invA;
   float CT[9];
   float detA;

   invA=(float *)calloc(9,sizeof(float));

   CT[0] = A[4]*A[8] - A[5]*A[7];
   CT[1] = A[2]*A[7] - A[1]*A[8];
   CT[2] = A[1]*A[5] - A[2]*A[4];
   CT[3] = A[5]*A[6] - A[3]*A[8];
   CT[4] = A[0]*A[8] - A[2]*A[6];
   CT[5] = A[2]*A[3] - A[0]*A[5];
   CT[6] = A[3]*A[7] - A[4]*A[6];
   CT[7] = A[1]*A[6] - A[0]*A[7];
   CT[8] = A[0]*A[4] - A[1]*A[3];

   detA=det3(A);

   if(detA!=0.0)
   {
      for(int i=0; i<9; i++)  
      invA[i]=CT[i]/detA;
   }

   return(invA);
}

void inv3(float *A, float *invA)
{
   float CT[9];
   float detA;

   CT[0] = A[4]*A[8] - A[5]*A[7];
   CT[1] = A[2]*A[7] - A[1]*A[8];
   CT[2] = A[1]*A[5] - A[2]*A[4];
   CT[3] = A[5]*A[6] - A[3]*A[8];
   CT[4] = A[0]*A[8] - A[2]*A[6];
   CT[5] = A[2]*A[3] - A[0]*A[5];
   CT[6] = A[3]*A[7] - A[4]*A[6];
   CT[7] = A[1]*A[6] - A[0]*A[7];
   CT[8] = A[0]*A[4] - A[1]*A[3];

   detA=det3(A);

   if(detA!=0.0)
   {
      for(int i=0; i<9; i++)  
         invA[i]=CT[i]/detA;
   }
}

// inverse is computed as the transpose of the cofactors divided by the determinant
float *inv2(float *A)
{
	float *invA;
	float detA;

	invA=(float *)calloc(4,sizeof(float));

	detA=A[0]*A[3]-A[1]*A[2];

	if(detA!=0.0)
	{
		invA[0] = A[3]/detA;
		invA[1] = -A[1]/detA;
		invA[2] = -A[2]/detA;
		invA[3] = A[0]/detA;
	}

   	return(invA);
}

double *inv2(double *A)
{
	double *invA;
	double detA;

	invA=(double *)calloc(4,sizeof(double));

	detA=A[0]*A[3]-A[1]*A[2];

	if(detA!=0.0)
	{
		invA[0] = A[3]/detA;
		invA[1] = -A[1]/detA;
		invA[2] = -A[2]/detA;
		invA[3] = A[0]/detA;
	}

   	return(invA);
}

double *inv4(double *A)
{
   int i;

   double detA;
   double *invA;
   double C[16];
   double B[9];
   double *transC;

   invA=(double *)calloc(16,sizeof(double));

   B[0]=A[5];  B[1]=A[6];  B[2]=A[7]; 
   B[3]=A[9];  B[4]=A[10]; B[5]=A[11]; 
   B[6]=A[13]; B[7]=A[14]; B[8]=A[15]; 
   C[0]=det3(B);

   B[0]=A[4];  B[1]=A[6];  B[2]=A[7]; 
   B[3]=A[8];  B[4]=A[10]; B[5]=A[11]; 
   B[6]=A[12]; B[7]=A[14]; B[8]=A[15]; 
   C[1]=-det3(B);

   B[0]=A[4];  B[1]=A[5];   B[2]=A[7]; 
   B[3]=A[8];  B[4]=A[9];  B[5]=A[11]; 
   B[6]=A[12]; B[7]=A[13]; B[8]=A[15]; 
   C[2]=det3(B);

   B[0]=A[4];  B[1]=A[5];   B[2]=A[6]; 
   B[3]=A[8];  B[4]=A[9];  B[5]=A[10]; 
   B[6]=A[12]; B[7]=A[13]; B[8]=A[14]; 
   C[3]=-det3(B);

   B[0]=A[1];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[9];  B[4]=A[10]; B[5]=A[11];
   B[6]=A[13]; B[7]=A[14]; B[8]=A[15];
   C[4]=-det3(B);

   B[0]=A[0];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[8];  B[4]=A[10]; B[5]=A[11];
   B[6]=A[12]; B[7]=A[14]; B[8]=A[15];
   C[5]=det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[3];
   B[3]=A[8];  B[4]=A[9]; B[5]=A[11];
   B[6]=A[12]; B[7]=A[13]; B[8]=A[15];
   C[6]=-det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[2];
   B[3]=A[8];  B[4]=A[9];  B[5]=A[10];
   B[6]=A[12]; B[7]=A[13]; B[8]=A[14];
   C[7]=det3(B);

   B[0]=A[1];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[5];  B[4]=A[6];  B[5]=A[7];
   B[6]=A[13];  B[7]=A[14];  B[8]=A[15];
   C[8]=det3(B);

   B[0]=A[0];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[4];  B[4]=A[6];  B[5]=A[7];
   B[6]=A[12];  B[7]=A[14];  B[8]=A[15];
   C[9]=-det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[3];
   B[3]=A[4];  B[4]=A[5];  B[5]=A[7];
   B[6]=A[12];  B[7]=A[13];  B[8]=A[15];
   C[10]=det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[2];
   B[3]=A[4];  B[4]=A[5];  B[5]=A[6];
   B[6]=A[12];  B[7]=A[13];  B[8]=A[14];
   C[11]=-det3(B);

   B[0]=A[1];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[5];  B[4]=A[6];  B[5]=A[7];
   B[6]=A[9]; B[7]=A[10]; B[8]=A[11];
   C[12]=-det3(B);

   B[0]=A[0];  B[1]=A[2];  B[2]=A[3];
   B[3]=A[4];  B[4]=A[6];  B[5]=A[7];
   B[6]=A[8]; B[7]=A[10]; B[8]=A[11];
   C[13]=det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[3];
   B[3]=A[4];  B[4]=A[5];  B[5]=A[7];
   B[6]=A[8]; B[7]=A[9]; B[8]=A[11];
   C[14]=-det3(B);

   B[0]=A[0];  B[1]=A[1];  B[2]=A[2];
   B[3]=A[4];  B[4]=A[5];  B[5]=A[6];
   B[6]=A[8]; B[7]=A[9]; B[8]=A[10];
   C[15]=det3(B);

   transC=trans(C,4,4);

   detA=det4(A);

   if(detA!=0.0)
   {
      for(i=0;i<16;i++)  
         invA[i]=transC[i]/detA;
   }

   free(transC);

   return(invA);
}

//////////////////////////////////////////////////////////////////

double euclideandistance(float *r0, float *r1, int n)
{
	double sum=0.0;

	for(int i=0; i<n; i++)
		sum += (r1[i]-r0[i])*(r1[i]-r0[i]);

	sum = sqrt(sum);

	return(sum);
}

//////////////////////////////////////////////////////////////////

void crossProduct(float *a, float *b, float *c)
{
	float dum[3];

	dum[0] = a[1]*b[2] - a[2]*b[1];
	dum[1] = a[2]*b[0] - a[0]*b[2];
	dum[2] = a[0]*b[1] - a[1]*b[0];

	c[0] = dum[0];
	c[1] = dum[1];
	c[2] = dum[2];
}

void crossProduct(double *a, double *b, double *c)
{
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];

   return;
}

void getcomplement (double *U, double *V, double *W)
{
   double norm;

   norm = normalize(U,3);

   if(norm < ESMALL) return;

   if ( fabs(U[0]) >= fabs(U[1]) )
   {
      double invLength = 1/sqrt(U[0]*U[0] + U[2]*U[2]);
      V[0] = -U[2]*invLength;
      V[1] = 0.0;
      V[2] = U[0]*invLength;
      W[0] = U[1]*V[2];
      W[1] = U[2]*V[0] - U[0]*V[2];
      W[2] = -U[1]*V[0];
   }
   else
   {
      double invLength = 1/sqrt(U[1]*U[1] + U[2]*U[2]);
      V[0] = 0.0;
      V[1] = U[2]*invLength;
      V[2] = -U[1]*invLength;
      W[0] = U[1]*V[2] - U[2]*V[1];
      W[1] = -U[0]*V[2];
      W[2] = U[0]*V[1];
   }

   return;
}

void getcomplement (float *U, float *V, float *W)
{
   float norm;

   norm = normalize(U,3);

   if(norm < ESMALL) return;

   if ( fabs(U[0]) >= fabsf(U[1]) )
   {
      float invLength = 1/sqrtf(U[0]*U[0] + U[2]*U[2]);
      V[0] = -U[2]*invLength;
      V[1] = 0.0;
      V[2] = U[0]*invLength;
      W[0] = U[1]*V[2];
      W[1] = U[2]*V[0] - U[0]*V[2];
      W[2] = -U[1]*V[0];
   }
   else
   {
      float invLength = 1/sqrtf(U[1]*U[1] + U[2]*U[2]);
      V[0] = 0.0;
      V[1] = U[2]*invLength;
      V[2] = -U[1]*invLength;
      W[0] = U[1]*V[2] - U[2]*V[1];
      W[1] = -U[0]*V[2];
      W[2] = U[0]*V[1];
   }

   return;
}

int ComputeRank(double *M)
{
   // Compute the maximum magnitude matrix entry.
   double abs, save, max = -1.0;
   int row, col, maxRow = -1, maxCol = -1;

   for(row = 0; row <3; row++)
   {
      for (col=row; col<3; col++)
      {
         abs = fabs( M[row*3 + col] );
         if (abs > max)
         {
            max = abs;
            maxRow = row;
            maxCol = col;
         }
      }
   }

   if (max < ESMALL)
   {
      // The rank is 0. The eigenvalue has multiplicity 3.
      return 0;
   }

   // The rank is at least 1. Swap the row containing the maximum-magnitude
   // entry with row 0.
   if (maxRow != 0)
   {
      for (col=0; col<3; col++)
      {
         save = M[col];
         M[col] = M[maxRow*3 + col];
         M[maxRow*3 + col] = save;
      }
   }

   // Row-reduce the matrix...

   // Scale row 0 to generate a 1-valued pivot.
   double invMax = 1/M[maxCol];
   M[0] *= invMax;
   M[1] *= invMax;
   M[2] *= invMax;

   // Eliminate the maxCol column entries in rows 1 and 2.
   if (maxCol == 0)
   {
      M[4] -= M[3]*M[1];
      M[5] -= M[3]*M[2];
      M[7] -= M[6]*M[1];
      M[8] -= M[6]*M[2];
      M[3] = 0.0;
      M[6] = 0.0;
   }
   else if (maxCol == 1)
   {
      M[3] -= M[4]*M[0];
      M[5] -= M[4]*M[2];
      M[6] -= M[7]*M[0];
      M[8] -= M[7]*M[2];
      M[4] = 0.0;
      M[7] = 0.0;
   }
   else
   {
      M[3] -= M[5]*M[0];
      M[4] -= M[5]*M[1];
      M[6] -= M[8]*M[0];
      M[7] -= M[8]*M[1];
      M[5] = 0.0;
      M[8] = 0.0;
   }

   // Compute the maximum-magnitude entry of the last two rows of the
   // row-reduced matrix.
   max = -1.0;
   maxRow = -1;
   maxCol = -1;
   for (row = 1; row < 3; row++)
   {
      for (col = 0; col < 3; col++)
      {
         abs = fabs( M[row*3 + col] );
         if (abs > max)
         {
            max = abs;
            maxRow = row;
            maxCol = col;
         }
      }
   }

   if (max < ESMALL)
   {
      // The rank is 1. The eigenvalue has multiplicity 2.
      return 1;
   }

   // If row 2 has the maximum-magnitude entry, swap it with row 1.
   if (maxRow == 2)
   {
      save = M[3];
      M[3] = M[6];
      M[6] = save;

      save = M[4];
      M[4] = M[7];
      M[7] = save;

      save = M[5];
      M[5] = M[8];
      M[8] = save;
   }

   // Scale row 1 to generate a 1-valued pivot.
   invMax = 1/M[3+maxCol];
   M[3] *= invMax;
   M[4] *= invMax;
   M[5] *= invMax;
   // The rank is 2. The eigenvalue has multiplicity 1.

   return 2;
}

int ComputeRank(float *M)
{
   // Compute the maximum magnitude matrix entry.
   float abs, save, max = -1.0;
   int row, col, maxRow = -1, maxCol = -1;

   for(row = 0; row <3; row++)
   {
      for (col=row; col<3; col++)
      {
         abs = fabsf( M[row*3 + col] );
         if (abs > max)
         {
            max = abs;
            maxRow = row;
            maxCol = col;
         }
      }
   }

   if (max < ESMALL)
   {
      // The rank is 0. The eigenvalue has multiplicity 3.
      return 0;
   }

   // The rank is at least 1. Swap the row containing the maximum-magnitude
   // entry with row 0.
   if (maxRow != 0)
   {
      for (col=0; col<3; col++)
      {
         save = M[col];
         M[col] = M[maxRow*3 + col];
         M[maxRow*3 + col] = save;
      }
   }

   // Row-reduce the matrix...

   // Scale row 0 to generate a 1-valued pivot.
   float invMax = 1/M[maxCol];
   M[0] *= invMax;
   M[1] *= invMax;
   M[2] *= invMax;

   // Eliminate the maxCol column entries in rows 1 and 2.
   if (maxCol == 0)
   {
      M[4] -= M[3]*M[1];
      M[5] -= M[3]*M[2];
      M[7] -= M[6]*M[1];
      M[8] -= M[6]*M[2];
      M[3] = 0.0;
      M[6] = 0.0;
   }
   else if (maxCol == 1)
   {
      M[3] -= M[4]*M[0];
      M[5] -= M[4]*M[2];
      M[6] -= M[7]*M[0];
      M[8] -= M[7]*M[2];
      M[4] = 0.0;
      M[7] = 0.0;
   }
   else
   {
      M[3] -= M[5]*M[0];
      M[4] -= M[5]*M[1];
      M[6] -= M[8]*M[0];
      M[7] -= M[8]*M[1];
      M[5] = 0.0;
      M[8] = 0.0;
   }

   // Compute the maximum-magnitude entry of the last two rows of the
   // row-reduced matrix.
   max = -1.0;
   maxRow = -1;
   maxCol = -1;
   for (row = 1; row < 3; row++)
   {
      for (col = 0; col < 3; col++)
      {
         abs = fabsf( M[row*3 + col] );
         if (abs > max)
         {
            max = abs;
            maxRow = row;
            maxCol = col;
         }
      }
   }

   if (max < ESMALL)
   {
      // The rank is 1. The eigenvalue has multiplicity 2.
      return 1;
   }

   // If row 2 has the maximum-magnitude entry, swap it with row 1.
   if (maxRow == 2)
   {
      save = M[3];
      M[3] = M[6];
      M[6] = save;

      save = M[4];
      M[4] = M[7];
      M[7] = save;

      save = M[5];
      M[5] = M[8];
      M[8] = save;
   }

   // Scale row 1 to generate a 1-valued pivot.
   invMax = 1/M[3+maxCol];
   M[3] *= invMax;
   M[4] *= invMax;
   M[5] *= invMax;
   // The rank is 2. The eigenvalue has multiplicity 1.

   return 2;
}

// A is a 3x3 symmetric matrix in 6x1 vector form
void s3eigenvec(double *A, double *evalue, double *UT)
{
   double M0[9], M1[9];
   int rank0, rank1;
   double row0[3], row1[3];
   double L;

   // M0 = A - evalue[0]*I; 
   L = evalue[0];
   M0[0] = A[0]-L;
   M0[4] = A[1]-L;
   M0[8] = A[2]-L;
   M0[3] = M0[1] = A[3];
   M0[7] = M0[5] = A[4];
   M0[6] = M0[2] = A[5];

   rank0 = ComputeRank(M0);
   if (rank0 == 0)
   {
      // evalue[0] = evalue[1] = evalue[2]
      for(int i=0; i<9; i++) UT[i] = 0.0;
      UT[0] = 1;
      UT[4] = 1;
      UT[8] = 1;

      return;
   }

   if (rank0 == 1)
   {
      // evalue[0] = evalue[1] > evalue[2]

      for(int i=0; i<3; i++) row0[i]=M0[i];

      getcomplement(row0,UT,UT+3);
      crossProduct(UT, UT+3, UT+6);

      return;
   }

   // rank0 == 2
   for(int i=0; i<3; i++) 
   {   
      row0[i]=M0[i];
      row1[i]=M0[3+i];
   }

   crossProduct(row0, row1, UT);
   normalize(UT,3);

   //M1 = A - evalue[1]*I;
   L = evalue[1];
   M1[0] = A[0]-L;
   M1[4] = A[1]-L;
   M1[8] = A[2]-L;
   M1[3] = M1[1] = A[3];
   M1[7] = M1[5] = A[4];
   M1[6] = M1[2] = A[5];

   rank1 = ComputeRank(M1);

   if (rank1 == 1)
   {
      // evalue[0] > evalue[1] = evalue[2]
      getcomplement(UT,UT+3,UT+6);
      return;
   }

   // rank1 == 2
   for(int i=0; i<3; i++) 
   {
      row0[i]=M1[i];
      row1[i]=M1[3+i];
   }

   crossProduct(row0, row1, UT+3);
   normalize(UT+3,3);

   crossProduct(UT, UT+3, UT+6);

   return;
}

void s3eigenvec(float *A, float *evalue, float *UT)
{
   float M0[9], M1[9];
   int rank0, rank1;
   float row0[3], row1[3];
   float L;

   // M0 = A - evalue[0]*I; 
   L = evalue[0];
   M0[0] = A[0]-L;
   M0[4] = A[1]-L;
   M0[8] = A[2]-L;
   M0[3] = M0[1] = A[3];
   M0[7] = M0[5] = A[4];
   M0[6] = M0[2] = A[5];

   rank0 = ComputeRank(M0);
   if (rank0 == 0)
   {
      // evalue[0] = evalue[1] = evalue[2]
      for(int i=0; i<9; i++) UT[i] = 0.0;
      UT[0] = 1;
      UT[4] = 1;
      UT[8] = 1;

      return;
   }

   if (rank0 == 1)
   {
      // evalue[0] = evalue[1] > evalue[2]

      for(int i=0; i<3; i++) row0[i]=M0[i];

      getcomplement(row0,UT,UT+3);
      crossProduct(UT, UT+3, UT+6);

      return;
   }

   // rank0 == 2
   for(int i=0; i<3; i++) 
   {   
      row0[i]=M0[i];
      row1[i]=M0[3+i];
   }

   crossProduct(row0, row1, UT);
   normalize(UT,3);

   //M1 = A - evalue[1]*I;
   L = evalue[1];
   M1[0] = A[0]-L;
   M1[4] = A[1]-L;
   M1[8] = A[2]-L;
   M1[3] = M1[1] = A[3];
   M1[7] = M1[5] = A[4];
   M1[6] = M1[2] = A[5];

   rank1 = ComputeRank(M1);

   if (rank1 == 1)
   {
      // evalue[0] > evalue[1] = evalue[2]
      getcomplement(UT,UT+3,UT+6);
      return;
   }

   // rank1 == 2
   for(int i=0; i<3; i++) 
   {
      row0[i]=M1[i];
      row1[i]=M1[3+i];
   }

   crossProduct(row0, row1, UT+3);
   normalize(UT+3,3);

   crossProduct(UT, UT+3, UT+6);

   return;
}

// L: 3x1 eigenvalues
// UT: 3x3 (row are eigenvectors)
// ULUT: U*L*U^{T} given as a 6x1 vector
void s3ULUT(double *L, double *UT, double *ULUT)
{
   ULUT[0] = UT[0]*UT[0]*L[0] + UT[3]*UT[3]*L[1] + UT[6]*UT[6]*L[2];
   ULUT[1] = UT[1]*UT[1]*L[0] + UT[4]*UT[4]*L[1] + UT[7]*UT[7]*L[2];
   ULUT[2] = UT[2]*UT[2]*L[0] + UT[5]*UT[5]*L[1] + UT[8]*UT[8]*L[2];

   ULUT[3] = UT[0]*UT[1]*L[0] + UT[3]*UT[4]*L[1] + UT[6]*UT[7]*L[2];
   ULUT[4] = UT[1]*UT[2]*L[0] + UT[4]*UT[5]*L[1] + UT[7]*UT[8]*L[2];
   ULUT[5] = UT[0]*UT[2]*L[0] + UT[3]*UT[5]*L[1] + UT[6]*UT[8]*L[2];
}

void s3ULUT(float *L, float *UT, float *ULUT)
{
   ULUT[0] = UT[0]*UT[0]*L[0] + UT[3]*UT[3]*L[1] + UT[6]*UT[6]*L[2];
   ULUT[1] = UT[1]*UT[1]*L[0] + UT[4]*UT[4]*L[1] + UT[7]*UT[7]*L[2];
   ULUT[2] = UT[2]*UT[2]*L[0] + UT[5]*UT[5]*L[1] + UT[8]*UT[8]*L[2];

   ULUT[3] = UT[0]*UT[1]*L[0] + UT[3]*UT[4]*L[1] + UT[6]*UT[7]*L[2];
   ULUT[4] = UT[1]*UT[2]*L[0] + UT[4]*UT[5]*L[1] + UT[7]*UT[8]*L[2];
   ULUT[5] = UT[0]*UT[2]*L[0] + UT[3]*UT[5]*L[1] + UT[6]*UT[8]*L[2];
}

// Computes A*B*A
// A, B, and ABA are 3x3 symmetric matrices represented by 6x1 vectors
void s3ABA(double *A, double *B, double *ABA)
{
   double ba00, ba01, ba02;
   double ba10, ba11, ba12;
   double ba20, ba21, ba22;

   double aba0, aba1, aba2, aba3, aba4, aba5;

   ba00 = B[0]*A[0] + B[3]*A[3] + B[5]*A[5];
   ba01 = B[0]*A[3] + B[3]*A[1] + B[5]*A[4];
   ba02 = B[0]*A[5] + B[3]*A[4] + B[5]*A[2];

   ba10 = B[3]*A[0] + B[1]*A[3] + B[4]*A[5];
   ba11 = B[3]*A[3] + B[1]*A[1] + B[4]*A[4];
   ba12 = B[3]*A[5] + B[1]*A[4] + B[4]*A[2];

   ba20 = B[5]*A[0] + B[4]*A[3] + B[2]*A[5];
   ba21 = B[5]*A[3] + B[4]*A[1] + B[2]*A[4];
   ba22 = B[5]*A[5] + B[4]*A[4] + B[2]*A[2];

   aba0 = A[0]*ba00 + A[3]*ba10 + A[5]*ba20;
   aba1 = A[3]*ba01 + A[1]*ba11 + A[4]*ba21;
   aba2 = A[5]*ba02 + A[4]*ba12 + A[2]*ba22;
   aba3 = A[3]*ba00 + A[1]*ba10 + A[4]*ba20;
   aba4 = A[5]*ba01 + A[4]*ba11 + A[2]*ba21;
   aba5 = A[5]*ba00 + A[4]*ba10 + A[2]*ba20;

   ABA[0] = aba0;
   ABA[1] = aba1;
   ABA[2] = aba2;
   ABA[3] = aba3;
   ABA[4] = aba4;
   ABA[5] = aba5;

   return;
}

void s3ABA(float *A, float *B, float *ABA)
{
   float ba00, ba01, ba02;
   float ba10, ba11, ba12;
   float ba20, ba21, ba22;

   float aba0, aba1, aba2, aba3, aba4, aba5;

   ba00 = B[0]*A[0] + B[3]*A[3] + B[5]*A[5];
   ba01 = B[0]*A[3] + B[3]*A[1] + B[5]*A[4];
   ba02 = B[0]*A[5] + B[3]*A[4] + B[5]*A[2];

   ba10 = B[3]*A[0] + B[1]*A[3] + B[4]*A[5];
   ba11 = B[3]*A[3] + B[1]*A[1] + B[4]*A[4];
   ba12 = B[3]*A[5] + B[1]*A[4] + B[4]*A[2];

   ba20 = B[5]*A[0] + B[4]*A[3] + B[2]*A[5];
   ba21 = B[5]*A[3] + B[4]*A[1] + B[2]*A[4];
   ba22 = B[5]*A[5] + B[4]*A[4] + B[2]*A[2];

   aba0 = A[0]*ba00 + A[3]*ba10 + A[5]*ba20;
   aba1 = A[3]*ba01 + A[1]*ba11 + A[4]*ba21;
   aba2 = A[5]*ba02 + A[4]*ba12 + A[2]*ba22;
   aba3 = A[3]*ba00 + A[1]*ba10 + A[4]*ba20;
   aba4 = A[5]*ba01 + A[4]*ba11 + A[2]*ba21;
   aba5 = A[5]*ba00 + A[4]*ba10 + A[2]*ba20;

   ABA[0] = aba0;
   ABA[1] = aba1;
   ABA[2] = aba2;
   ABA[3] = aba3;
   ABA[4] = aba4;
   ABA[5] = aba5;

   return;
}

// Computes the square of the Riemammian distance between matrices D and F
// given D and invsqrtF = F^{-1/2}
double p3RiemannianDistance(double *D, double *invsqrtF)
{
   double M[6]; // M = F^{-1/2} * D * F^{-1/2}
   double L[3];

   s3ABA(invsqrtF, D, M);
   s3eigenval(M,L);

   assert(L[0]>ESMALL && L[1]>ESMALL && L[2]>ESMALL);

   return( log(L[0])*log(L[0]) + log(L[1])*log(L[1]) + log(L[2])*log(L[2]) );
}

float p3RiemannianDistance(float *D, float *invsqrtF)
{
   float M[6]; // M = F^{-1/2} * D * F^{-1/2}
   float L[3];

   s3ABA(invsqrtF, D, M);
   s3eigenval(M,L);

   if( L[0]<ESMALL || L[1]<ESMALL || L[2]<ESMALL )
      return(-1.0);

   return( logf(L[0])*logf(L[0]) + logf(L[1])*logf(L[1]) + logf(L[2])*logf(L[2]) );
}

// L are the eigenvalues of D*F^{-1}
double p3RiemannianDistance(double *L)
{
   assert(L[0]>ESMALL && L[1]>ESMALL && L[2]>ESMALL);

   return( log(L[0])*log(L[0]) + log(L[1])*log(L[1]) + log(L[2])*log(L[2]) );
}

// A and invsqrtA are 3x3 symmetric positive definite matrices
// represented as 6x1 vectors
void p3invsqrt(double *A, double *invsqrtA)
{
   double L[3];
   double UT[9];

   // compute eigenvalues of A
   s3eigenval(A,L);

   // compute eigenvectors of A
   s3eigenvec(A,L,UT);

   assert(L[0]>ESMALL && L[1]>ESMALL && L[2]>ESMALL);

   // compute eigenvalues of A^{-1/2}
   for(int i=0; i<3; i++) L[i] = 1.0/sqrt(L[i]);

   // compute A^{-1/2}
   s3ULUT(L, UT, invsqrtA);

   return;
}

// A and invsqrtA are 3x3 symmetric positive definite matrices
// represented as 6x1 vectors
int p3invsqrt(float *A, float *invsqrtA)
{
   float L[3];
   float UT[9];

   // compute eigenvalues of A
   s3eigenval(A,L);

   // compute eigenvectors of A
   s3eigenvec(A,L,UT);

   if( L[0]<ESMALL || L[1]<ESMALL || L[2]<ESMALL )
      return 0;

   // compute eigenvalues of A^{-1/2}
   for(int i=0; i<3; i++) L[i] = 1.0/sqrtf(L[i]);

   // compute A^{-1/2}
   s3ULUT(L, UT, invsqrtA);

   return 1;
}

void p3invsqrt(double *A, double *invsqrtA, double *sqrtA)
{
   double L[3];
   double UT[9];

   // compute eigenvalues of A
   s3eigenval(A,L);

   // compute eigenvectors of A
   s3eigenvec(A,L,UT);

//   assert(L[0]>ESMALL && L[1]>ESMALL && L[2]>ESMALL);

   // compute eigenvalues of A^{1/2}
   for(int i=0; i<3; i++) L[i] = sqrt(L[i]);

   // compute A^{1/2}
   s3ULUT(L, UT, sqrtA);

   // compute eigenvalues of A^{1/2}
   for(int i=0; i<3; i++) L[i] = 1.0/L[i];

   // compute A^{-1/2}
   s3ULUT(L, UT, invsqrtA);

   return;
}

// WRONG! Under construction
void s3multi(double *A, double *B, double *AB)
{
   AB[0] = B[0]*A[0] + B[3]*A[3] + B[5]*A[5];
   AB[1] = B[0]*A[3] + B[3]*A[1] + B[5]*A[4];
   AB[2] = B[0]*A[5] + B[3]*A[4] + B[5]*A[2];

   AB[3] = B[3]*A[0] + B[1]*A[3] + B[4]*A[5];
   AB[4] = B[3]*A[3] + B[1]*A[1] + B[4]*A[4];
   AB[5] = B[3]*A[5] + B[1]*A[4] + B[4]*A[2];

   AB[6] = B[5]*A[0] + B[4]*A[3] + B[2]*A[5];
   AB[7] = B[5]*A[3] + B[4]*A[1] + B[2]*A[4];
   AB[8] = B[5]*A[5] + B[4]*A[4] + B[2]*A[2];

   return;
}

void p3update(double *D, double *invsqrtD, double *sqrtD, double *W, double epsilon)
{
   double E[6];
   double L[3];
   double UT[9];
   double EXPW[6];

   s3ABA(invsqrtD, W, E);
   s3eigenval(E, L);
   s3eigenvec(E, L, UT);
   for(int i=0; i<3; i++) L[i] = exp(-epsilon*L[i]);
   s3ULUT(L, UT, EXPW);
   s3ABA(sqrtD, EXPW, D);
}

void p3update(double *D, double *W, double epsilon)
{
   double invsqrtD[6], sqrtD[6];
   double E[6];
   double L[3];
   double UT[9];
   double EXPW[6];

   p3invsqrt(D, invsqrtD, sqrtD);
   s3ABA(invsqrtD, W, E);
   s3eigenval(E, L);
   s3eigenvec(E, L, UT);
   for(int i=0; i<3; i++) L[i] = exp(-epsilon*L[i]);
   s3ULUT(L, UT, EXPW);
   s3ABA(sqrtD, EXPW, D);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// qqT: Input 6x1 vector representing the 3x3 symmetrix matrix q*transpose(q).
// D: Input 6x1 vector representing 3x3 symmetrix matrix D.
// Output: Scalar value of matrix multiplication transpose(q)*D*q.
// Note: The output is the Frobenius inner product of 3x3 matrices q*transpose(q) and D.
//////////////////////////////////////////////////////////////////////////////////////////////////
double Frobenius_s3(double *qqT, float *D)
{
   return( 
      qqT[0]*D[0] + 
      qqT[1]*D[1] + 
      qqT[2]*D[2] + 
      qqT[3]*D[3]*2.0 + 
      qqT[4]*D[4]*2.0 + 
      qqT[5]*D[5]*2.0 
   );
}

double Frobenius_s3(double *qqT, double *D)
{
   return( 
      qqT[0]*D[0] + 
      qqT[1]*D[1] + 
      qqT[2]*D[2] + 
      qqT[3]*D[3]*2.0 + 
      qqT[4]*D[4]*2.0 + 
      qqT[5]*D[5]*2.0 
   );
}
