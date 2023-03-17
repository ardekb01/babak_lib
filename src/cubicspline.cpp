#define _cubicspline

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

///////////////////////////////////
// For Cubic B-spline interpolation
///////////////////////////////////
#define K0 10
#define MAXINDEX 100000
///////////////////////////////////

float *computeBeta(float *del)
{
	float *beta;
	float x,dum;

	*del = 2.0/MAXINDEX;
	beta = (float *)calloc(MAXINDEX+1,sizeof(float));
	for(int i=0; i<=MAXINDEX; i++)
	{
		x = (*del) * i;
		if( x<=1.0 )
		{
			dum = x*x;
			beta[i] = dum*x*3.0 - 6.0*dum + 4.0;
		}
		else
		{
			dum = -x + 2.0;
			beta[i] = dum*dum*dum;
		}
	}

	return(beta);	
}

float cubicSplineSynthesis(float *c, int nx, int ny, int nz, float x, float y, float z, float *beta, float del)
{
	double ux,uy,uz;
	float betax[4], betay[4], betaz[4]; 
	float Ldel;

	float sum;
	int Kx, Ky, Kz, J;
	int np;
	int slice_ptr, row_ptr;

	np = nx*ny;

	if(x<0.0 || y<0.0 || z<0.0 || x>(nx-1.0) || y>(ny-1) || z>(nz-1) )
		return(0.0);
	
	ux = (x - (int)x + 1.0)/del;
	uy = (y - (int)y + 1.0)/del;
	uz = (z - (int)z + 1.0)/del;

	for(int L=0; L<=3; L++)
	{
		Ldel = L/del;

		J = (int)rint( ux - Ldel ); 	
		if (J<0) J = -J;
		if (J>MAXINDEX) J = MAXINDEX;
		betax[L] = beta[J];

		J = (int)rint( uy - Ldel ); 	
		if (J<0) J = -J;
		if (J>MAXINDEX) J = MAXINDEX;
		betay[L] = beta[J];

		J = (int)rint( uz - Ldel ); 	
		if (J<0) J = -J;
		if (J>MAXINDEX) J = MAXINDEX;
		betaz[L] = beta[J];
	}

	sum = 0.0;
	for(int Lz=0; Lz<=3; Lz++)
	{
		Kz = Lz + (int)z - 1;

		// ensure Kz is in the range from 0 to nz-1
		if(Kz<0) Kz = -Kz; else if(Kz>=nz) Kz=2*nz-2-Kz;

		slice_ptr = Kz*np;

		for(int Ly=0; Ly<=3; Ly++)
		{
			Ky = Ly + (int)y - 1;
			
			// ensure Ky is in the range from 0 to ny-1
			if(Ky<0) Ky = -Ky; else if(Ky>=ny) Ky=2*ny-2-Ky;

			row_ptr = slice_ptr + Ky*nx;

			for(int Lx=0; Lx<=3; Lx++)
			{
				Kx = Lx + (int)x - 1;

				// ensure Kx is in the range from 0 to nx-1
				if(Kx<0) Kx = -Kx; else if(Kx>=nx) Kx=2*nx-2-Kx;

				sum += c[row_ptr + Kx]*betax[Lx]*betay[Ly]*betaz[Lz];
			}
		}
	}

	return(sum);
}

// Implements the algorithm in Box 2
// The coefficients are returned in s.
void cubicSplineAnalysis(float *s, int N)
{
	float z1;
	float *cplus;
	float *cminus;
	float sum;
	float z1_k; 		// z1 to the power of k

	z1 = (float)(sqrt(3.0) - 2.0);

	cplus = s;

	// compute c+(0)
	sum=0.0;
	z1_k=1.0; 
	for(int k=0; k<=K0; k++)
	{
		sum += (z1_k * s[k]);
		z1_k *= z1;
	}
	cplus[0] = sum;

	// compute c+(k) k=1, ..., N-1
	for(int k=1; k<N; k++)
		cplus[k] = s[k] + z1*cplus[k-1];

	cminus = cplus;

	// compute c-(N-1)
	cminus[N-1] = z1*(cplus[N-1] + z1*cplus[N-2])/(z1*z1 - 1.0);

	// compute c-(k) k=N-2, N-3, ..., 0
	for(int k=(N-2); k>=0; k--)
		cminus[k] = z1 * (cminus[k+1] - cplus[k] );

	// The factor of 6 is deliberately left out
	// for(int k=0; k<N; k++) cminus[k] *= 6.0;

	return;
}

void cubicSplineAnalysis(unsigned char *s, float *c, int N)
{
	float z1;
	float *cplus;
	float *cminus;
	float sum;
	float z1_k; 		// z1 to the power of k

	z1 = (float)(sqrt(3.0) - 2.0);

	cplus = c;

	// compute c+(0)
	sum=0.0;
	z1_k=1.0; 
	for(int k=0; k<=K0; k++)
	{
		sum += (z1_k * s[k]);
		z1_k *= z1;
	}
	cplus[0] = sum;

	// compute c+(k) k=1, ..., N-1
	for(int k=1; k<N; k++)
		cplus[k] = s[k] + z1*cplus[k-1];

	cminus = cplus;

	// compute c-(N-1)
	cminus[N-1] = z1*(cplus[N-1] + z1*cplus[N-2])/(z1*z1 - 1.0);

	// compute c-(k) k=N-2, N-3, ..., 0
	for(int k=(N-2); k>=0; k--)
		cminus[k] = z1 * (cminus[k+1] - cplus[k] );

	// The factor of 6 is deliberately left out
	// for(int k=0; k<N; k++) cminus[k] *= 6.0;

	return;
}

void cubicSplineAnalysis(short *s, float *c, int N)
{
	float z1;
	float *cplus;
	float *cminus;
	float sum;
	float z1_k; 		// z1 to the power of k

	z1 = (float)(sqrt(3.0) - 2.0);

	cplus = c;

	// compute c+(0)
	sum=0.0;
	z1_k=1.0; 
	for(int k=0; k<=K0; k++)
	{
		sum += (z1_k * s[k]);
		z1_k *= z1;
	}
	cplus[0] = sum;

	// compute c+(k) k=1, ..., N-1
	for(int k=1; k<N; k++)
		cplus[k] = s[k] + z1*cplus[k-1];

	cminus = cplus;

	// compute c-(N-1)
	cminus[N-1] = z1*(cplus[N-1] + z1*cplus[N-2])/(z1*z1 - 1.0);

	// compute c-(k) k=N-2, N-3, ..., 0
	for(int k=(N-2); k>=0; k--)
		cminus[k] = z1 * (cminus[k+1] - cplus[k] );

	// The factor of 6 is deliberately left out
	// for(int k=0; k<N; k++) cminus[k] *= 6.0;

	return;
}

void cubicSplineAnalysis(float *s, float *c, int N)
{
	float z1;
	float *cplus;
	float *cminus;
	float sum;
	float z1_k; 		// z1 to the power of k

	z1 = (float)(sqrt(3.0) - 2.0);

	cplus = c;

	// compute c+(0)
	sum=0.0;
	z1_k=1.0; 
	for(int k=0; k<=K0; k++)
	{
		sum += (z1_k * s[k]);
		z1_k *= z1;
	}
	cplus[0] = sum;

	// compute c+(k) k=1, ..., N-1
	for(int k=1; k<N; k++)
		cplus[k] = s[k] + z1*cplus[k-1];

	cminus = cplus;

	// compute c-(N-1)
	cminus[N-1] = z1*(cplus[N-1] + z1*cplus[N-2])/(z1*z1 - 1.0);

	// compute c-(k) k=N-2, N-3, ..., 0
	for(int k=(N-2); k>=0; k--)
		cminus[k] = z1 * (cminus[k+1] - cplus[k] );

	// The factor of 6 is deliberately left out
	// for(int k=0; k<N; k++) cminus[k] *= 6.0;

	return;
}

void cubicSplineAnalysis(unsigned char *s, float *c, int nx, int ny)
{
	float *dum;
	int ptr;

	dum = (float *)calloc(ny,sizeof(float));

	// do every row
	for(int j=0; j<ny; j++)
	{
		ptr = j*nx;

		cubicSplineAnalysis(s+ptr, c+ptr, nx);
	}

	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++) dum[j] = c[j*nx+ i];
		cubicSplineAnalysis(dum, ny);
		for(int j=0; j<ny; j++) c[j*nx + i] = dum[j];
	}

	free(dum);
}

void cubicSplineAnalysis(short *s, float *c, int nx, int ny)
{
	float *dum;
	int ptr;

	dum = (float *)calloc(ny,sizeof(float));

	// do every row
	for(int j=0; j<ny; j++)
	{
		ptr = j*nx;

		cubicSplineAnalysis(s+ptr, c+ptr, nx);
	}

	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++) dum[j] = c[j*nx+ i];
		cubicSplineAnalysis(dum, ny);
		for(int j=0; j<ny; j++) c[j*nx + i] = dum[j];
	}

	free(dum);
}

void cubicSplineAnalysis(float *s, float *c, int nx, int ny)
{
	float *dum;
	int ptr;

	dum = (float *)calloc(ny,sizeof(float));

	// do every row
	for(int j=0; j<ny; j++)
	{
		ptr = j*nx;

		cubicSplineAnalysis(s+ptr, c+ptr, nx);
	}

	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++) dum[j] = c[j*nx+ i];
		cubicSplineAnalysis(dum, ny);
		for(int j=0; j<ny; j++) c[j*nx + i] = dum[j];
	}

	free(dum);
}

void cubicSplineAnalysis(unsigned char *s, float *c, int nx, int ny, int nz)
{
	float *dum;
	int np;
	int column_ptr, row_ptr, slice_ptr;

	np = nx*ny;

	dum = (float *)calloc(nz,sizeof(float));

	// do every slice
	for(int k=0; k<nz; k++)
	{
		slice_ptr = k*np;
		cubicSplineAnalysis(s+slice_ptr, c+slice_ptr, nx, ny);
	}

	for(int j=0; j<ny; j++)
	{
		row_ptr = j*nx;

		for(int i=0; i<nx; i++)
		{
			column_ptr = row_ptr + i;

			for(int k=0; k<nz; k++) dum[k] = c[k*np + column_ptr];

			cubicSplineAnalysis(dum, nz);

			for(int k=0; k<nz; k++) c[k*np + column_ptr] = dum[k];
		}
	}

	free(dum);
}

void cubicSplineAnalysis(short *s, float *c, int nx, int ny, int nz)
{
	float *dum;
	int np;
	int column_ptr, row_ptr, slice_ptr;

	np = nx*ny;

	dum = (float *)calloc(nz,sizeof(float));

	// do every slice
	for(int k=0; k<nz; k++)
	{
		slice_ptr = k*np;
		cubicSplineAnalysis(s+slice_ptr, c+slice_ptr, nx, ny);
	}

	for(int j=0; j<ny; j++)
	{
		row_ptr = j*nx;

		for(int i=0; i<nx; i++)
		{
			column_ptr = row_ptr + i;

			for(int k=0; k<nz; k++) dum[k] = c[k*np + column_ptr];

			cubicSplineAnalysis(dum, nz);

			for(int k=0; k<nz; k++) c[k*np + column_ptr] = dum[k];
		}
	}

	free(dum);
}

void cubicSplineAnalysis(float *s, float *c, int nx, int ny, int nz)
{
	float *dum;
	int np;
	int column_ptr, row_ptr, slice_ptr;

	np = nx*ny;

	dum = (float *)calloc(nz,sizeof(float));

	// do every slice
	for(int k=0; k<nz; k++)
	{
		slice_ptr = k*np;
		cubicSplineAnalysis(s+slice_ptr, c+slice_ptr, nx, ny);
	}

	for(int j=0; j<ny; j++)
	{
		row_ptr = j*nx;

		for(int i=0; i<nx; i++)
		{
			column_ptr = row_ptr + i;

			for(int k=0; k<nz; k++) dum[k] = c[k*np + column_ptr];

			cubicSplineAnalysis(dum, nz);

			for(int k=0; k<nz; k++) c[k*np + column_ptr] = dum[k];
		}
	}

	free(dum);
}

short *resliceImageCubicSpline(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1, 
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T)
{
	float *beta, del;
	float *c;
	float xc1,yc1,zc1; 
	float xc2,yc2,zc2; 
	float Ax,Bx;
	float Ay,By;
	float Az,Bz;
  	float xx,yy,zz; /* translation parameters */
  	float x,y,z;   

	int i,j,k;     /* loop index */

	int q;
	int np2,nv2;

	short *im2;

	beta = computeBeta(&del);
	
    c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
    cubicSplineAnalysis(im1, c, nx1, ny1, nz1);

//////////////////////////////////////////////////////

	np2=nx2*ny2;
	nv2=np2*nz2;

	im2=(short *)calloc(nv2,sizeof(short));

	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

	xc1=dx1*(nx1-1)/2.0;      /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

	T[0] /= dx1;
	T[1] /= dx1;
	T[2] /= dx1;
	T[3] /= dx1;
	T[3] += xc1/dx1;

	T[4] /= dy1;
	T[5] /= dy1;
	T[6] /= dy1;
	T[7] /= dy1;
	T[7] += yc1/dy1;

	T[8]  /= dz1;
	T[9]  /= dz1;
	T[10] /= dz1;
	T[11] /= dz1;
	T[11] += zc1/dz1;

	q=0;
	for(k=0;k<nz2;k++) 
	{
  		zz=k*dz2-zc2;
	  	Bx=T[2]*zz+T[3];
	  	By=T[6]*zz+T[7];
		Bz=T[10]*zz+T[11];
		for(j=0;j<ny2;j++) 
		{
     		yy=j*dy2-yc2;
			Ax=T[1]*yy+Bx;
			Ay=T[5]*yy+By;
			Az=T[9]*yy+Bz;

  			for(i=0;i<nx2;i++) 
			{
        		xx=i*dx2-xc2;

           		x=T[0]*xx+Ax;
	    	   	y=T[4]*xx+Ay;
	    	   	z=T[8]*xx+Az;

				im2[q++]= (short)(cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
			}
		}
	}

	free(beta);
	free(c);
	return( im2 );
}

float *resliceImageCubicSpline(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1, 
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T)
{
	float *beta, del;
	float *c;
	float xc1,yc1,zc1; 
	float xc2,yc2,zc2; 
	float Ax,Bx;
	float Ay,By;
	float Az,Bz;
  	float xx,yy,zz; /* translation parameters */
  	float x,y,z;   

	int i,j,k;     /* loop index */

	int q;
	int np2,nv2;

	float *im2;

	beta = computeBeta(&del);
	
	c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
	cubicSplineAnalysis(im1, c, nx1, ny1, nz1);

//////////////////////////////////////////////////////

	np2=nx2*ny2;
	nv2=np2*nz2;

	im2=(float *)calloc(nv2,sizeof(float));

	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

	xc1=dx1*(nx1-1)/2.0;      /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

	T[0] /= dx1;
	T[1] /= dx1;
	T[2] /= dx1;
	T[3] /= dx1;
	T[3] += xc1/dx1;

	T[4] /= dy1;
	T[5] /= dy1;
	T[6] /= dy1;
	T[7] /= dy1;
	T[7] += yc1/dy1;

	T[8]  /= dz1;
	T[9]  /= dz1;
	T[10] /= dz1;
	T[11] /= dz1;
	T[11] += zc1/dz1;

	q=0;
	for(k=0;k<nz2;k++) 
	{
  		zz=k*dz2-zc2;
	  	Bx=T[2]*zz+T[3];
	  	By=T[6]*zz+T[7];
		Bz=T[10]*zz+T[11];
		for(j=0;j<ny2;j++) 
		{
     			yy=j*dy2-yc2;
			Ax=T[1]*yy+Bx;
			Ay=T[5]*yy+By;
			Az=T[9]*yy+Bz;

  			for(i=0;i<nx2;i++) 
			{
        		xx=i*dx2-xc2;

           		x=T[0]*xx+Ax;
	    	   	y=T[4]*xx+Ay;
	    	   	z=T[8]*xx+Az;

				im2[q++]= (float)(cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
			}
		}
	}

	free(beta);
	free(c);
	return( im2 );
}
