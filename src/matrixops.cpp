#define _matrixops

#include <stdio.h>
#include <stdlib.h>

#include <f2c.h>
#include <clapack.h>

#include "babak_lib.h"
#include "bbk_linear_algebra.h"

void solveATxEqB(float *A, int nr, int nc, float *B);

float *ATB_float(float *a,int nra,int nca,float *b,int nrb,int ncb);

/* diagATA_float:
   Diagnolized the (n x n) symmetric matrix ATA.
   n is the order of ATA matrix
   ATA is n*n block of memory
   Set uplo='L' or 'l' when the input matrix ATA lower triangular,
   and uplo='U' or 'u' when the input matrix ATA upper triangular.
   Returns the eigenvalues (diagonal elemenets) in
   descending order.   Contents of ATA are replaced
   by U' where ATA = U D U' and D is the diagonal
   matrix.
*/
float *diagATA_float(float *ATA, int n, char uplo);

/* ATA_float:
   Computes and returns the symmetric (nc X nc) matrix A'A.
   nr is the number of rows in A
   nc is the number of columns of A
   A is nr*nc block of memory
   If uplo='L' or 'l' the result is lower triangular,
   otherwise the result is upper triangular.
*/
float *ATA_float(float *A,int nr,int nc, char uplo);

/* AAT_float:
   Computes and returns the symmetric matrix (nr X nr) AA'.
   nr is the number of rows in A
   nc is the number of columns of A
   A is nr*nc block of memory
   If uplo='L' or 'l' the result is lower triangular,
   otherwise the result is upper triangular.
*/
float *AAT_float(float *A,int nr,int nc, char uplo);

/* Ax_float:
   Computes and returns the vector (nr X 1) Ax.
   nr is the number of rows in A
   nc is the number of columns of A
   A is nr*nc block of memory
   x is nc*1 block of memory
*/
float *Ax_float(float *A,int nr,int nc, float *x);

/* ATx_float:
   Computes and returns the vector (nc X 1) A'x.
   nr is the number of rows in A
   nc is the number of columns of A
   A is nr*nc block of memory
   x is nr*1 block of memory
*/
float *ATx_float(float *A,int nr,int nc, float *x);

float *ATA_float(float *A,int nr,int nc, char uplo)
{
   char trans='N';

   float alpha=1.0;
   float beta=0.0;
   float *ATA;

   ATA=(float *)calloc(nc*nc,sizeof(float));

   if(ATA==NULL) return(NULL);

   if(uplo=='L' || uplo=='l' ) uplo='U';
   else (uplo='L');

   ssyrk_(&uplo, &trans, (integer *)&nc, (integer *)&nr, &alpha, A, (integer *)&nc, &beta, ATA, (integer *)&nc);

   return(ATA);
}

float *AAT_float(float *A,int nr,int nc, char uplo)
{
   char trans='T';

   float alpha=1.0;
   float beta=0.0;
   float *AAT;

   AAT=(float *)calloc(nr*nr,sizeof(float));

   if(AAT==NULL) return(NULL);

   if(uplo=='L' || uplo=='l' ) uplo='U';
   else (uplo='L');

   ssyrk_(&uplo,&trans,(integer *)&nr,(integer *)&nc,&alpha,A,(integer *)&nc,&beta,AAT,(integer *)&nr);

   return(AAT);
}

float *Ax_float(float *A,int nr,int nc, float *x)
{
   char trans='T';

   float *Ax;
   float alpha=1.0;
   float beta=0.0;

   int inc=1;

   Ax=(float *)calloc(nr,sizeof(float));
   if(Ax==NULL) return(NULL);

   sgemv_(&trans,(integer *)&nc,(integer *)&nr,&alpha,A,(integer *)&nc,x,(integer *)&inc,&beta,Ax,(integer *)&inc);

   return(Ax);
}

float *ATx_float(float *A,int nr,int nc, float *x)
{
   char trans='N';

   float *ATx;
   float alpha=1.0;
   float beta=0.0;

   int inc=1;

   ATx=(float *)calloc(nc,sizeof(float));
   if(ATx==NULL) return(NULL);

   sgemv_(&trans,(integer *)&nc,(integer *)&nr,&alpha,A,(integer *)&nc,x,(integer *)&inc,&beta,ATx,(integer *)&inc);

   return(ATx);
}

float *diagATA_float(float *ATA, int n, char uplo)
{
   int i,j;
   char JOBZ='V';
   float *D;
   float *WORK;
   float dum;
   float *dummy;
   int LWORK;
   int INFO;
   int m;
   
   if(uplo=='L' || uplo=='l') uplo='U';
   else uplo='L';

   D=(float *)calloc(n,sizeof(float));

   LWORK=3*n;
   WORK=(float *)calloc(LWORK,sizeof(float));

   ssyev_(&JOBZ,&uplo,(integer *)&n, ATA, (integer *)&n, D, WORK, (integer *)&LWORK,(integer *)&INFO );

   free(WORK);

   /* reverse order of elements of D */
   m=(int)(n/2);
   for(i=0;i<m;i++)
   {
      dum=D[i];   
      D[i]=D[n-1-i];
      D[n-1-i]=dum;
   }

   /* reverse order of rows of ATA */
   dummy=(float *)calloc(n,sizeof(float));
   for(i=0;i<m;i++)
   {
      for(j=0;j<n;j++)
         dummy[j]=ATA[i*n + j];   

      for(j=0;j<n;j++)
         ATA[i*n + j]=ATA[(n-1-i)*n + j];   

      for(j=0;j<n;j++)
         ATA[(n-1-i)*n + j]=dummy[j];   
   }

   free(dummy);

   if(INFO==0)
      return(D);
   else 
      return(NULL);
}

float *ATB_float(float *A,int ,int ncA,float *B,int nrB,int ncB)
{
   char trans1='N';
   char trans2='T';

   float *arr1,*arr2;

   float *c;
   float alpha=1.0;
   float beta=0.0;

   int nr1,nc1,nr2;

   arr1=B;
   arr2=A;

   nr1=ncB; nc1=nrB;
   nr2=ncA;

   c=(float *)calloc(nr1*nr2,sizeof(float));
   if(c==NULL) return(NULL);

   sgemm_(&trans1,&trans2,(integer *)&nr1,(integer *)&nr2,(integer *)&nc1,&alpha,arr1,(integer *)&nr1,arr2,(integer *)&nr2,&beta,c,(integer *)&nr1);

   return(c);
}

void solveATxEqB(float *A, int nr, int nc, float *B)
{
   int INFO;
   int NRHS=1;
   int *IPIV;

   IPIV=(int *)calloc(nr,sizeof(int));

   sgesv_((integer *)&nr, (integer *)&NRHS, A, (integer *)&nc, (integer *)IPIV, B, (integer *)&nc, (integer *)&INFO );

   free(IPIV);
}

///////////////////////////////////////////////////////////////////////////
void mat_mat_trans(float *A,int Ar,int Ac,float *B,int Br, float *C)
{
   float sum;
 
   if(Ar <= 0 || Ac <= 0 || Br <= 0)
   {
      fprintf(stderr,"mat_mat_trans(): Negative or zero dimensions!\n");
      return;
   }
 
   for(int i=0;i<Ar;i++)
   for(int j=0;j<Br;j++)
   {
      sum=0.0;
      for(int k=0;k<Ac;k++)
         sum += A[i*Ac+k]*B[j*Ac+k];
 
      C[i*Br+j]=sum;
   }
 
   return;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void mat_trans_mat(float *A, int Ar, int Ac, float *B, int Bc, float *C)
{
	float sum;
 
	if(Ar <= 0 || Ac <= 0 || Bc <= 0)
	{
		fprintf(stderr,"mat_trans_mat(): Negative or zero dimensions!\n");
		return;
	}
 
	for(int i=0; i<Ac; i++)
	for(int j=0; j<Bc; j++)
	{
		sum=0.0;
		for(int k=0;k<Ar;k++)
			sum += A[k*Ac+i]*B[k*Bc+j];
 
		C[i*Bc+j]=sum;
	}
 
	return;
}

void mat_trans_mat(double *A, int Ar, int Ac, double *B, int Bc, double *C)
{
	double sum;
 
	if(Ar <= 0 || Ac <= 0 || Bc <= 0)
	{
		fprintf(stderr,"mat_trans_mat(): Negative or zero dimensions!\n");
		return;
	}
 
	for(int i=0; i<Ac; i++)
	for(int j=0; j<Bc; j++)
	{
		sum=0.0;
		for(int k=0;k<Ar;k++)
			sum += A[k*Ac+i]*B[k*Bc+j];
 
		C[i*Bc+j]=sum;
	}
 
	return;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

// Computes At = V * S * Ut
// A = MxN matrix (M>N). At=transpose(A) is an NxM matrix.
// M = Number of rows of A (Number of columns of At)
// N = Number of columns of A (Number of rows of At)
// Ut = NxM matrix = transpose(U)
// S = N vector
// V = NxN matrix
// IMPORTANT: contents of At will be destroyed.
void svd(float *At, int m, int n, float *Ut, float *V, float *S)
// void svd(float *At, integer M, integer N, float *Ut, float *V, float *S)
{
	integer M;
	integer N;
	char JOBU='S';
	char JOBVT='A';

	integer LWORK=2048;
	integer INFO;

	float WORK[2048];

	M = m;
	N = n;

	sgesvd_(&JOBU,&JOBVT,&M,&N,At,&M, S, Ut,&M,V,&N,WORK,&LWORK,&INFO );
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////

float *projectionMatrix(float *X, int N, int p, int *rank)
{
	int r;
	float *Ut, *V, *S;
	float *Xt;
	float *P;

	// memory allocation
	P = (float *)calloc(N*N,sizeof(float));

	// memory allocation
	Xt = (float *)calloc(p*N,sizeof(float));

	for(int i=0; i<N*p; i++) Xt[i]=X[i];
	transpose_matrix(Xt, N,  p);

	// memory allocation
	Ut = (float *)calloc(p*N,sizeof(float));

	// memory allocation
	V = (float *)calloc(p*p,sizeof(float));

	// memory allocation
	S = (float *)calloc(p,sizeof(float));

	svd(Xt, N, p, Ut, V, S);

    // Determine rank of X based on it's number of 'non-zero' singular values
    r=0;
    for(int i=0; i<p; i++)
    {
        if( (int)(1000.0*S[i]/S[0]) != 0)
            r++;
        else
            S[i]=0.0;
    }

	// Set the last (p-r) rows of Ut to 0.
	for(int i=r; i<p; i++)
	for(int j=0; j<N; j++) Ut[i*N + j]=0.0;

	mat_trans_mat(Ut, p, N, Ut, N, P);

	*rank = r;

	free(Xt);
	free(S);
	free(Ut);
	free(V);

	return(P);
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

float *projectionMatrix(double *X, int N, int p, int *rank)
{
	int r;
	float *Ut, *V, *S;
	float *Xt;
	float *P;

	// memory allocation
	P = (float *)calloc(N*N,sizeof(float));

	// memory allocation
	Xt = (float *)calloc(p*N,sizeof(float));

	for(int i=0; i<N*p; i++) Xt[i]=(float)X[i];
	transpose_matrix(Xt, N,  p);

	// memory allocation
	Ut = (float *)calloc(p*N,sizeof(float));

	// memory allocation
	V = (float *)calloc(p*p,sizeof(float));

	// memory allocation
	S = (float *)calloc(p,sizeof(float));

	svd(Xt, N, p, Ut, V, S);

    // Determine rank of X based on it's number of 'non-zero' singular values
    r=0;
    for(int i=0; i<p; i++)
    {
        if( (int)(1000.0*S[i]/S[0]) != 0)
            r++;
        else
            S[i]=0.0;
    }

	// Set the last (p-r) rows of Ut to 0.
	for(int i=r; i<p; i++)
	for(int j=0; j<N; j++) Ut[i*N + j]=0.0;

	mat_trans_mat(Ut, p, N, Ut, N, P);

	*rank = r;

	free(Xt);
	free(S);
	free(Ut);
	free(V);

	return(P);
}

void projectVector(double *x, double *xpar, double *xper, float *Pz, int n)
{
	multi(Pz, n, n, x, n, 1, xpar);
	for(int i=0; i<n; i++) xper[i] = x[i] - xpar[i];
}

// Input A is a square NxN matrix
// N specifies the size of matrix A
// n is between 0 and N-1 and specifies the row and column index to be zeroed.
int zeroRowCol(float *A, int N, int n)
{
	if(A==NULL) return(1);
	if(N<=0) return(2);
	if(n>=N || n<0) return(3);

	for(int j=0; j<N; j++) 
		A[n*N + j] = 0.0;
		
	for(int i=0; i<N; i++) 
		A[i*N + n] = 0.0;

	return(0);
}

// Input A is a square NxN matrix
// N specifies the size of matrix A
// n is between 0 and N-1 and specifies the row and column index to be set.
// a is a vector of size N which replaces the nth row and nth column of A.
int setRowCol(float *A, int N, int n, float *a)
{
	if(A==NULL) return(1);
	if(N<=0) return(2);
	if(n>=N || n<0) return(3);
	if(a==NULL) return(4);

	for(int j=0; j<N; j++) 
		A[n*N + j] = a[j];
		
	for(int i=0; i<N; i++) 
		A[i*N + n] = a[i];

	return(0);
}

void zeroVector(float *v, int n)
{
	for(int i=0; i<n; i++) v[i]=0.0;
}

void zeroVector(char *v, int n)
{
	for(int i=0; i<n; i++) v[i]=0;
}

void oneVector(float *v, int n)
{
	for(int i=0; i<n; i++) v[i]=1.0;
}

void oneVector(char *v, int n)
{
	for(int i=0; i<n; i++) v[i]=1;
}
