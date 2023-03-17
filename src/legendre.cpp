#define _legendre

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

double *LegendreAnalysis(float *image, int nx, int ny, int nz, int mx, int my, int mz);
float *LegendreSynthesis(double *c, int nx, int ny, int nz, int mx, int my, int mz);

void LegendreSynthesis(double *c, int mx, int my, float *image, int nx, int ny);

void LegendrePoly(double *p0, double *q0, double *p1, double *q1, double *x, int N, int n);
void LegendrePoly(float *p0, float *q0, float *p1, float *q1, float *x, int N, int n);
void integral_1d(double *a, short *b, int n, double *d);
void integral_1d(double *a, float *b, int n, double *d);
void integral_1d(double *a, double *b, int n, double *d);

////////////////////////////////////////////////////////////////////////////////////
// Returns an (n+1)xN matrix. The ith row of this matrix contains the Legendre
// polynomial of degree i-1. Each polynomial is evaluated at N uniformly spaced
// points from -1.0 to 1.0 (including points -1.0 and 1.0). These points are
// given by the second row of the matrix since P1(x)=x.
////////////////////////////////////////////////////////////////////////////////////
double *LegendrePolynomials(int n, int N)
{
	double *P;
	double *x;
	double *p_1, *p_2, *p;

	// requires N to be greater than 1
	if(N<=1) return(NULL);

	// requires n to be greater than or equal to 0
	if(n<0) return(NULL);

	x = (double *)calloc(N, sizeof(double));
	if(x==NULL) return(NULL);

	P = (double *)calloc( (n+1)*N ,sizeof(double));
	if(P==NULL) return(NULL);

	// x is an array of N real numbers uniformly spaced from 
	// -1.0 to 1.0 with x[0]=-1.0 and x[N-1]=1.0;
	for(int i=0; i<N; i++) x[i] = -1.0 + i*2.0/(N-1.0);

	if(n==0) 	// form P0 
	{
		for(int i=0; i<N; i++) P[i] = 1.0;
	}
	else if(n==1) 	// form P0 and P1 
	{
		for(int i=0; i<N; i++) 
		{
			P[i] = 1.0;
			P[i+N] = x[i];
		}
	}
	else
	{
		// form P0 and P1
		for(int i=0; i<N; i++) 
		{
			P[i] = 1.0;
			P[i+N] = x[i];
		}

		// form P2, ..., Pn
		for(int j=2; j<=n; j++) 
		{
			p = P + j*N;			// p points to the current row of matrix P
			p_1 = P + (j-1)*N;		// p_1 points to the row above p
			p_2 = P + (j-2)*N;		// p_2 points to the row above p_1

			for(int i=0; i<N; i++) 
				p[i] = ( (2*j-1)*x[i]*p_1[i] - (j-1)*p_2[i] )/j;
		}
	}

	free(x);

	return(P);
}
////////////////////////////////////////////////////////////////////////////////////
// end of LegendrePolynomials()
////////////////////////////////////////////////////////////////////////////////////

void integral_1d(double *a, short *b, int n, double *d)
{
	double sum=0.0;

	for(int i=0; i<n; i++)
		sum += a[i]*b[i];

	*d = sum*2.0/n;

	return;
}

void integral_1d(double *a, float *b, int n, double *d)
{
	double sum=0.0;

	for(int i=0; i<n; i++)
		sum += a[i]*b[i];

	*d = sum*2.0/n;

	return;
}

void integral_1d(double *a, double *b, int n, double *d)
{
	double sum=0.0;

	for(int i=0; i<n; i++)
		sum += a[i]*b[i];

	*d = sum*2.0/n;

	return;
}

// Inputs:
// p0 : Legendre polynomical of order n-1 (evaluated at the N points in x)
// q0 : first derivative p0
// n : order of the output Legendre polynomial
// N : dimension of the arrays p0,q0,p1,q1, and x.
// x : abscissas
// Outputs:
// p1 : Legendre polynomical of order n 
// q1 : derivative of p1

void LegendrePoly(double *p0, double *q0, double *p1, double *q1, double *x, int N, int n)
{
	for(int i=0; i<N; i++)
	{
		q1[i] = x[i]*q0[i] + n*p0[i];
		p1[i] = x[i]*p0[i] + (x[i]*x[i]-1.0)*q0[i]/n;
	}

	return;
}

void LegendrePoly(float *p0, float *q0, float *p1, float *q1, float *x, int N, int n)
{
	for(int i=0; i<N; i++)
	{
		q1[i] = x[i]*q0[i] + n*p0[i];
		p1[i] = x[i]*p0[i] + (x[i]*x[i]-1.0)*q0[i]/n;
	}

	return;
}

#if 0
int main(int argc, char **argv)
{
	FILE *fp1,*fp2;
	int np;

	float dx,dy,dz;

	int mx=16;
	int nx;

	int my=16;
	int ny;

	int mz=16;
	int nz;

	double *c;			// series coefficients
	float *image;
	float *image_est;

	fp1=fopen(argv[1],"r");
	fp2=fopen("tt.wrp","w");

	fread(&nx,sizeof(int),1,fp1);
	fread(&ny,sizeof(int),1,fp1);
	fread(&nz,sizeof(int),1,fp1);
	fread(&dx,sizeof(float),1,fp1);
	fread(&dy,sizeof(float),1,fp1);
	fread(&dz,sizeof(float),1,fp1);

	printf("%d %d %d\n",nx,ny,nz);
	printf("%f %f %f\n",dx,dy,dz);

	fwrite(&nx,sizeof(int),1,fp2);
	fwrite(&ny,sizeof(int),1,fp2);
	fwrite(&nz,sizeof(int),1,fp2);
	fwrite(&dx,sizeof(float),1,fp2);
	fwrite(&dy,sizeof(float),1,fp2);
	fwrite(&dz,sizeof(float),1,fp2);

	image = (float *)calloc(nx*ny*nz,sizeof(float));

	fread(image,sizeof(float),nx*ny*nz,fp1);
	c=LegendreAnalysis(image, nx, ny, nz, mx, my, mz);
	image_est=LegendreSynthesis(c, nx, ny, nz, mx, my, mz);
	fwrite(image_est,sizeof(float),nx*ny*nz,fp2);
	free(c); free(image_est);
	printf("OK1\n");

	fread(image,sizeof(float),nx*ny*nz,fp1);
	c=LegendreAnalysis(image, nx, ny, nz, mx, my, mz);
	image_est=LegendreSynthesis(c, nx, ny, nz, mx, my, mz);
	fwrite(image_est,sizeof(float),nx*ny*nz,fp2);
	free(c); free(image_est);
	printf("OK2\n");

	fread(image,sizeof(float),nx*ny*nz,fp1);
	c=LegendreAnalysis(image, nx, ny, nz, mx, my, mz);
	image_est=LegendreSynthesis(c, nx, ny, nz, mx, my, mz);
	fwrite(image_est,sizeof(float),nx*ny*nz,fp2);
	free(c); free(image_est);
	printf("OK3\n");

	fclose(fp1);
	fclose(fp2);
	////////////////////////////////////////////////////////////////////////
	np = nx*ny;

/*
    fp=fopen("fit.txt","w");
    for(int i=0; i<nx; i++)
        fprintf(fp,"%d %lf %lf\n", i, image[91*np+55*nx+i] ,image_est[91*np+55*nx+i]);
	fclose(fp);
*/

/*
	fp=fopen("xw_est.img","w");
	for(int i=0; i<nx*ny*nz; i++)
		temp[i] = (short)(image_est[i]+0.5);
	fwrite(temp,sizeof(short),nx*ny*nz,fp);
	fclose(fp);
*/
}
#endif

double *LegendreAnalysis(float *image, int nx, int ny, int nz, int mx, int my, int mz)
{
   double *x;
   double *Px,*Qx;
   double *y;
   double *Py,*Qy;
   double *z;
   double *Pz,*Qz;

   double *c;			// series coefficients
   int mp;
   int np;

   np = nx*ny;
   mp = mx*my;

   /////////////////////////////////////////////////////////////////////////
   c = (double *)calloc(mx*my*mz,sizeof(double));
   /////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////////
   x = (double *)calloc(nx,sizeof(double));
   Px = (double *)calloc(mx*nx,sizeof(double));
   Qx = (double *)calloc(mx*nx,sizeof(double));

   // x is an array of nx real number going from -(nx-1.0)/nx to (nx-1.0)/nx
   for(int i=0; i<nx; i++)
      x[i] = -(nx-1.0)/nx + i*2.0/nx;

   // initialization
   for(int i=0; i<nx; i++)
   {
      Px[i]=1.0; 
      Qx[i]=0.0;
   }

   for(int q=1; q<mx; q++)
      LegendrePoly(Px+(q-1)*nx, Qx+(q-1)*nx, Px+q*nx, Qx+q*nx, x, nx, q);

   for(int q=0; q<mx; q++)
   for(int i=0; i<nx; i++)
      Px[q*nx + i] /= sqrt( 2. / (2.0*q+1) );
   /////////////////////////////////////////////////////////////////////////

   y = (double *)calloc(ny,sizeof(double));
   Py = (double *)calloc(my*ny,sizeof(double));
   Qy = (double *)calloc(my*ny,sizeof(double));

   // set y
   for(int i=0; i<ny; i++)
   {
      y[i] = -(2.0/ny)*((ny-1)/2.0) + i*2.0/ny;
   }

   // initialization
   for(int i=0; i<ny; i++)
   {
      Py[i]=1.0; 
      Qy[i]=0.0;
   }

   for(int r=1; r<my; r++)
      LegendrePoly(Py+(r-1)*ny, Qy+(r-1)*ny, Py+r*ny, Qy+r*ny, y, ny, r);

   for(int r=0; r<my; r++)
   for(int i=0; i<ny; i++)
      Py[r*ny + i] /= sqrt( 2. / (2.0*r+1) );
   /////////////////////////////////////////////////////////////////////////

   z = (double *)calloc(nz,sizeof(double));
   Pz = (double *)calloc(mz*nz,sizeof(double));
   Qz = (double *)calloc(mz*nz,sizeof(double));

   // set z
   for(int i=0; i<nz; i++)
   {
      z[i] = -(2.0/nz)*((nz-1)/2.0) + i*2.0/nz;
   }

   // initialization
   for(int i=0; i<nz; i++)
   {
      Pz[i]=1.0; 
      Qz[i]=0.0;
   }

   for(int s=1; s<mz; s++)
      LegendrePoly(Pz+(s-1)*nz, Qz+(s-1)*nz, Pz+s*nz, Qz+s*nz, z, nz, s);

   for(int s=0; s<mz; s++)
   for(int i=0; i<nz; i++)
      Pz[s*nz + i] /= sqrt( 2. / (2.0*s+1) );
   /////////////////////////////////////////////////////////////////////////

   double *dum1d;
   double *dum2d;
   
   dum1d = (double *)calloc(nz,sizeof(double));
   dum2d = (double *)calloc(nz*ny,sizeof(double));

   for(int q=0; q<mx; q++)
   {
      for(int j=0; j<ny; j++)
      for(int k=0; k<nz; k++)
         integral_1d(Px+q*nx, image+k*np+j*nx, nx, dum2d+k*ny+j);

      for(int r=0; r<my; r++)
      {
         for(int k=0; k<nz; k++)
            integral_1d(Py+r*ny, dum2d+k*ny, ny, dum1d+k);

         for(int s=0; s<mz; s++)
            integral_1d(Pz+s*nz, dum1d, nz, c+s*mp+r*mx+q);
      }
   }
   free(dum1d);
   free(dum2d);
   ////////////////////////////////////////////////////////////////////////

   free(x); free(y); free(z);
   free(Px); free(Qx);
   free(Py); free(Qy);
   free(Pz); free(Qz);

   return(c);
}

float *LegendreSynthesis(double *c, int nx, int ny, int nz, int mx, int my, int mz)
{
	double *x;
	double *Px,*Qx;
	double *y;
	double *Py,*Qy;
	double *z;
	double *Pz,*Qz;

	float *image;
	double *sum1;
	double *sum2;

	int mp;
	int np;

	mp = mx*my;
	np = nx*ny;

	/////////////////////////////////////////////////////////////////////////
	x = (double *)calloc(nx,sizeof(double));
	Px = (double *)calloc(mx*nx,sizeof(double));
	Qx = (double *)calloc(mx*nx,sizeof(double));

	// set x
	for(int i=0; i<nx; i++)
	{
		x[i] = -(2.0/nx)*((nx-1)/2.0) + i*2.0/nx;
	}

	// initialization
	for(int i=0; i<nx; i++)
	{
		Px[i]=1.0; 
		Qx[i]=0.0;
	}

	for(int q=1; q<mx; q++)
		LegendrePoly(Px+(q-1)*nx, Qx+(q-1)*nx, Px+q*nx, Qx+q*nx, x, nx, q);

	for(int q=0; q<mx; q++)
	for(int i=0; i<nx; i++)
		Px[q*nx + i] /= sqrt( 2. / (2.0*q+1) );
	/////////////////////////////////////////////////////////////////////////

	y = (double *)calloc(ny,sizeof(double));
	Py = (double *)calloc(my*ny,sizeof(double));
	Qy = (double *)calloc(my*ny,sizeof(double));

	// set y
	for(int i=0; i<ny; i++)
	{
		y[i] = -(2.0/ny)*((ny-1)/2.0) + i*2.0/ny;
	}

	// initialization
	for(int i=0; i<ny; i++)
	{
		Py[i]=1.0; 
		Qy[i]=0.0;
	}

	for(int r=1; r<my; r++)
		LegendrePoly(Py+(r-1)*ny, Qy+(r-1)*ny, Py+r*ny, Qy+r*ny, y, ny, r);

	for(int r=0; r<my; r++)
	for(int i=0; i<ny; i++)
		Py[r*ny + i] /= sqrt( 2. / (2.0*r+1) );
	/////////////////////////////////////////////////////////////////////////

	z = (double *)calloc(nz,sizeof(double));
	Pz = (double *)calloc(mz*nz,sizeof(double));
	Qz = (double *)calloc(mz*nz,sizeof(double));

	// set z
	for(int i=0; i<nz; i++)
	{
		z[i] = -(2.0/nz)*((nz-1)/2.0) + i*2.0/nz;
	}

	// initialization
	for(int i=0; i<nz; i++)
	{
		Pz[i]=1.0; 
		Qz[i]=0.0;
	}

	for(int s=1; s<mz; s++)
		LegendrePoly(Pz+(s-1)*nz, Qz+(s-1)*nz, Pz+s*nz, Qz+s*nz, z, nz, s);

	for(int s=0; s<mz; s++)
	for(int i=0; i<nz; i++)
		Pz[s*nz + i] /= sqrt( 2. / (2.0*s+1) );
	/////////////////////////////////////////////////////////////////////////

	image = (float *)calloc(nx*ny*nz,sizeof(float));

	sum1 = (double *)calloc(mx*my*nz,sizeof(double));	

	for(int q=0; q<mx; q++)
	for(int r=0; r<my; r++)
	for(int z=0; z<nz; z++)
	for(int s=0; s<mz; s++)
		sum1[z*mp + r*mx + q] += c[s*mp+r*mx+q]*Pz[s*nz + z];

	sum2 = (double *)calloc(mx*ny*nz,sizeof(double));	

	for(int q=0; q<mx; q++)
	for(int z=0; z<nz; z++)
	for(int y=0; y<ny; y++)
	for(int r=0; r<my; r++)
		sum2[z*mx*ny + y*mx + q] += sum1[z*mp + r*mx + q]*Py[r*ny + y];

	free(sum1);

	for(int z=0; z<nz; z++)
	for(int y=0; y<ny; y++)
	for(int x=0; x<nx; x++)
	for(int q=0; q<mx; q++)
		image[z*np + y*nx + x] += sum2[z*mx*ny + y*mx + q]*Px[q*nx + x];

	free(sum2);

	free(x); free(y); free(z);
	free(Px); free(Qx);
	free(Py); free(Qy);
	free(Pz); free(Qz);

	return(image);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


double *LegendreAnalysis(float *image, int nx, int ny, int mx, int my)
{
   double *x;
   double *Px,*Qx;
   double *y;
   double *Py,*Qy;

   double *c;			// series coefficients

   /////////////////////////////////////////////////////////////////////////
   c = (double *)calloc(mx*my,sizeof(double));
   /////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////////
   x = (double *)calloc(nx,sizeof(double));
   Px = (double *)calloc(mx*nx,sizeof(double));
   Qx = (double *)calloc(mx*nx,sizeof(double));

   // x is an array of nx real number going from -(nx-1.0)/nx to (nx-1.0)/nx
   for(int i=0; i<nx; i++)
      x[i] = -(nx-1.0)/nx + i*2.0/nx;

   // initialization
   for(int i=0; i<nx; i++)
   {
      Px[i]=1.0; 
      Qx[i]=0.0;
   }

   for(int q=1; q<mx; q++)
      LegendrePoly(Px+(q-1)*nx, Qx+(q-1)*nx, Px+q*nx, Qx+q*nx, x, nx, q);

   for(int q=0; q<mx; q++)
   for(int i=0; i<nx; i++)
      Px[q*nx + i] /= sqrt( 2. / (2.0*q+1) );
   /////////////////////////////////////////////////////////////////////////

   y = (double *)calloc(ny,sizeof(double));
   Py = (double *)calloc(my*ny,sizeof(double));
   Qy = (double *)calloc(my*ny,sizeof(double));

   // set y
   for(int i=0; i<ny; i++)
      y[i] = -(ny-1.0)/ny + i*2.0/ny;

   // initialization
   for(int i=0; i<ny; i++)
   {
      Py[i]=1.0; 
      Qy[i]=0.0;
   }

   for(int r=1; r<my; r++)
      LegendrePoly(Py+(r-1)*ny, Qy+(r-1)*ny, Py+r*ny, Qy+r*ny, y, ny, r);

   for(int r=0; r<my; r++)
   for(int i=0; i<ny; i++)
      Py[r*ny + i] /= sqrt( 2. / (2.0*r+1) );
   /////////////////////////////////////////////////////////////////////////

   double *dum;

   dum = (double *)calloc(ny,sizeof(double));
   for(int q=0; q<mx; q++)
   {
      for(int j=0; j<ny; j++)
         integral_1d(Px+q*nx, image+j*nx, nx, dum+j);

      for(int r=0; r<my; r++)
         integral_1d(Py+r*ny, dum, ny, c+r*mx+q);
   }

   free(dum);
   ////////////////////////////////////////////////////////////////////////

   free(x); free(y);
   free(Px); free(Qx);
   free(Py); free(Qy);

   return(c);
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

void LegendreSynthesis(double *c, int mx, int my, float *image, int nx, int ny)
{
	double *x;
	double *Px,*Qx;
	double *y;
	double *Py,*Qy;

	double *sum;

	/////////////////////////////////////////////////////////////////////////
	x = (double *)calloc(nx,sizeof(double));
	Px = (double *)calloc(mx*nx,sizeof(double));
	Qx = (double *)calloc(mx*nx,sizeof(double));

	// set x
	for(int i=0; i<nx; i++)
	{
		x[i] = -(2.0/nx)*((nx-1)/2.0) + i*2.0/nx;
	}

	// initialization
	for(int i=0; i<nx; i++)
	{
		Px[i]=1.0; 
		Qx[i]=0.0;
	}

	for(int q=1; q<mx; q++)
		LegendrePoly(Px+(q-1)*nx, Qx+(q-1)*nx, Px+q*nx, Qx+q*nx, x, nx, q);

	for(int q=0; q<mx; q++)
	for(int i=0; i<nx; i++)
		Px[q*nx + i] /= sqrt( 2. / (2.0*q+1) );
	/////////////////////////////////////////////////////////////////////////

	y = (double *)calloc(ny,sizeof(double));
	Py = (double *)calloc(my*ny,sizeof(double));
	Qy = (double *)calloc(my*ny,sizeof(double));

	// set y
	for(int i=0; i<ny; i++)
	{
		y[i] = -(2.0/ny)*((ny-1)/2.0) + i*2.0/ny;
	}

	// initialization
	for(int i=0; i<ny; i++)
	{
		Py[i]=1.0; 
		Qy[i]=0.0;
	}

	for(int r=1; r<my; r++)
		LegendrePoly(Py+(r-1)*ny, Qy+(r-1)*ny, Py+r*ny, Qy+r*ny, y, ny, r);

	for(int r=0; r<my; r++)
	for(int i=0; i<ny; i++)
		Py[r*ny + i] /= sqrt( 2. / (2.0*r+1) );
	/////////////////////////////////////////////////////////////////////////

	sum = (double *)calloc(mx*ny,sizeof(double));	

	for(int q=0; q<mx; q++)
	for(int y=0; y<ny; y++)
	for(int r=0; r<my; r++)
		sum[y*mx + q] += c[r*mx + q]*Py[r*ny + y];

	for(int y=0; y<ny; y++)
	for(int x=0; x<nx; x++)
	for(int q=0; q<mx; q++)
		image[y*nx + x] += sum[y*mx + q]*Px[q*nx + x];

	free(sum);

	free(x); free(y);
	free(Px); free(Qx);
	free(Py); free(Qy);
}
