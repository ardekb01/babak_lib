///////////////////////////////////////////////////////////////// 
// brainwash.cxx 
// Copyright 2016 by Babak A. Ardekani 
// ALL RIGHTS RESERVED. No part of this program may be 
// used, transferred, or modified by any means without 
// prior written permission. Making copies of any part
// of this program for any purpose is a violation of
// copyright laws.
///////////////////////////////////////////////////////////////// 

#include <stdio.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <dirent.h>
#include <volume.h>
#include <spm_analyze.h>
#include <babak_lib.h>
#include <sph.h>
#include <smooth.h>
#include <minmax.h>
#include <interpolator.h>
#include <stats.h>
#include <string.h>

#define YES 1
#define NO 0
#define MAXITER 10
#define MAXR 15
#define MAXNATLAS 999

#define XMATRIXSIZE 255
#define YMATRIXSIZE 255
#define ZMATRIXSIZE 189
#define VOXELSIZE 1.0

//////////////////////////////////////////////////////////////////////////////////////////////////

float FWHM=1.0;
char atlasfilename[MAXNATLAS][64];

// multi-resolution image dimensions
DIM dim1, dim2, dim4, dim8;

int opt;

static struct option options[] =
{
   {"-version", 0, 'V'},
   {"-Version", 0, 'V'},

   {"-n", 1, 'n'},

   {"-iter8", 1, '8'},
   {"-iter4", 1, '4'},
   {"-iter2", 1, '2'},
   {"-iter1", 1, '1'},

   {"-r", 1, 'r'},
   {"-R", 1, 'R'},

   {"-thresh", 1, 't'},
   {"-threshold", 1, 't'},
   {"-t", 1, 't'},

   {"-i", 1, 'i'},
   {"-lm", 1, 'M'},

   {"-atlaslist", 1, 'l'},
   {"-list", 1, 'l'},
   {"-l", 1, 'l'},

   {"-atlasdir", 1, 'd'},
   {"-dir", 1, 'd'},
   {"-d", 1, 'd'},

   {"-ppm", 0, 'p'},
   {"-v", 0, 'v'},
   {"-verbose", 0, 'v'},
   {"-h",0,'h'},
   {"-help",0,'h'},
   {0, 0, 0}
};

int opt_cubicspline=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

extern float detect_lm(SPH &searchsph, SPH &testsph, SHORTIM testim, int lmcm[], SPH &refsph, int lm[]);

void print_help_and_exit()
{
   printf("\n\nUsage:\n"
   "\tbrainwash [-v or -verbose] [-h or -help] [-version]\n"
   "\t[-r <patch radius>] [-R <search radius>] [-t <percent threshold>]\n"
   "\t-i <image>\n\n" 

   "Required arguments:\n\n"
   "\t-i <image>\n"
   "\t\tSpecifies the image to be skull-stripped. <image> is expected\n"
   "\t\tto be of type short in NIFTI format.\n\n"

   "Optional arguments:\n\n"
   "\t-v or -verbose\n"
   "\t\tEnables verbose mode.\n\n"

   "\t-r <patch radius>\n"
   "\t\tSpecifies the patch radius (default=3 pixels).\n\n"

   "\t-R <search radius>\n"
   "\t\tSpecifies the search radius (default=3 pixels).\n\n"

   "\t-h or -help\n"
   "\t\tPrints help message.\n\n"

   "\t-version\n"
   "\t\tPrints program vesion.\n\n"
   );
   exit(0);
}

short *computeReslicedImage2(short *im1, DIM dim1, DIM dim2, float *Xwarp, float *Ywarp, float *Zwarp)
{
 	float  x,y,z;   
	int q;
	short *im2;
	float xc1, yc1, zc1;
	float xc2, yc2, zc2;

	float *beta=NULL, del;
	float *c=NULL;

	if(opt_cubicspline)
	{
		beta=computeBeta(&del);
		c = (float *)calloc(dim1.nv, sizeof(float));
		cubicSplineAnalysis(im1, c, dim1.nx, dim1.ny, dim1.nz);
	}

	im2=(short *)calloc(dim2.nv,sizeof(short));

	xc1=dim1.dx*(dim1.nx-1)/2.0;     /* +---+---+ */
	yc1=dim1.dy*(dim1.ny-1)/2.0;
	zc1=dim1.dz*(dim1.nz-1)/2.0;

   	xc2=dim2.dx*(dim2.nx-1)/2.0;     /* +---+---+ */
	yc2=dim2.dy*(dim2.ny-1)/2.0;
	zc2=dim2.dz*(dim2.nz-1)/2.0;

	q=0;
	for(int k=0;k<dim2.nz;k++) 
	for(int j=0;j<dim2.ny;j++) 
  	for(int i=0;i<dim2.nx;i++) 
	{
		z = (k*dim2.dz - zc2 + Zwarp[q] + zc1) /dim1.dz;
		y = (j*dim2.dy - yc2 + Ywarp[q] + yc1) /dim1.dy;
		x = (i*dim2.dx - xc2 + Xwarp[q] + xc1) /dim1.dx;

		if(opt_cubicspline)
			im2[q++] = (short)(cubicSplineSynthesis(c, dim1.nx, dim1.ny, dim1.nz, x, y, z, beta, del)+0.5);
		else
			im2[q++]=(short)(linearInterpolator(x,y,z,im1,dim1.nx,dim1.ny,dim1.nz,dim1.np)+0.5);
	}

	if(opt_cubicspline)
	{
		free(beta);
		free(c);
	}

	return( im2 );
}

void brainwashnlReg(short *im1, short *im2, short *msk, DIM dim, int r, int R, float *Xwarp, float *Ywarp, float *Zwarp, int S)
{
   float I[16];
   float *Xw, *Yw, *Zw;

   SPH im1sph(r);
   SPH im2sph(r);
   SPH searchsph(R);

   SHORTIM sim1;
   set_dim(sim1,dim);
   sim1.v=im1;

   SHORTIM sim2;
   set_dim(sim2,dim);
   sim2.v=im2;

   int v; // voxel index going from 0 to dim.nv-1
   int P[3], Q[3];

   Xw = (float *)calloc(dim.nv, sizeof(float));
   Yw = (float *)calloc(dim.nv, sizeof(float));
   Zw = (float *)calloc(dim.nv, sizeof(float));

   for(int i=0; i<dim.nv; i++) if(msk[i]==0) im1[i]=im2[i]=0;
    
   for(int k=0; k<dim.nz; k++)
   for(int j=0; j<dim.ny; j++)
   for(int i=0; i<dim.nx; i++)
   {
      v = k*dim.np + j*dim.nx + i;

      if(msk[v]==100 || msk[v]<=0 || im2[v]<=0)
      {
         Xw[v] = Yw[v] = Zw[v] = 0.0;
         continue;
      }

      im2sph.set(sim2,i,j,k);
      standardize(im2sph.v,im2sph.v,im2sph.n);

      P[0]=i; P[1]=j; P[2]=k;
      detect_lm(searchsph, im1sph, sim1, P, im2sph, Q);

      // make sure Q is within the mask
      if(msk[ Q[2]*dim.np + Q[1]*dim.nx + Q[0] ] != 0)
      {
         Xw[v] = (Q[0]-P[0])*dim.dx;
         Yw[v] = (Q[1]-P[1])*dim.dy;
         Zw[v] = (Q[2]-P[2])*dim.dz;
      }
   }

   float *Xww, *Yww, *Zww;
   Xww = (float *)calloc(dim.nv, sizeof(float));
   Yww = (float *)calloc(dim.nv, sizeof(float));
   Zww = (float *)calloc(dim.nv, sizeof(float));

   {
      SPH xsph(S);
      SPH ysph(S);
      SPH zsph(S);

      for(int k=0; k<dim.nz; k++)
      for(int j=0; j<dim.ny; j++)
      for(int i=0; i<dim.nx; i++)
      {
         v = k*dim.np + j*dim.nx + i;

         if(msk[v]==100 || msk[v]<=0 || im2[v]<=0)
         {
            Xww[v] = Yww[v] = Zww[v] = 0.0;
            continue;
         }

         xsph.set(Xw,dim.nx,dim.ny,dim.nz,i,j,k);
         ysph.set(Yw,dim.nx,dim.ny,dim.nz,i,j,k);
         zsph.set(Zw,dim.nx,dim.ny,dim.nz,i,j,k);

         Xww[v] = xsph.mean();
         Yww[v] = ysph.mean();
         Zww[v] = zsph.mean();
      }

      // replace Xw with Xww etc.
      free(Xw); free(Yw); free(Zw);
      Xw=Xww; Yw=Yww; Zw=Zww;
   }

   float *tmp;

   set_to_I(I,4); 
   tmp = resliceImage(Xw, dim, dim1, I, LIN); 
   free(Xw); Xw=tmp;

   set_to_I(I,4); 
   tmp = resliceImage(Yw, dim, dim1, I, LIN); 
   free(Yw); Yw=tmp;

   set_to_I(I,4); 
   tmp = resliceImage(Zw, dim, dim1, I, LIN); 
   free(Zw); Zw=tmp;

   for(int i=0; i<dim1.nv; i++) {
      Xwarp[i] += Xw[i];
      Ywarp[i] += Yw[i];
      Zwarp[i] += Zw[i];
   }

   free(Xw); free(Yw); free(Zw);
}

// finds an affine registration in "A" to takes points from im1 to im2 space
void affineReg(short *im1, short *im2, short *msk, DIM dim, int r, int R, float *A)
{
   SPH im1sph(r);
   SPH im2sph(r);
   SPH searchsph(R);

   SHORTIM sim1; // structured im1
   set_dim(sim1,dim);
   sim1.v=im1;

   SHORTIM sim2;
   set_dim(sim2,dim);
   sim2.v=im2;

   int *PT, *QT;
   char *flg;

   float *P, *Q;

   flg = (char *)calloc(dim.nv, sizeof(char));
   PT = (int *)calloc(dim.nv*3, sizeof(int));
   QT = (int *)calloc(dim.nv*3, sizeof(int));
   
   int c=0;
   int v;
   for(int k=0; k<dim.nz; k++)
   for(int j=0; j<dim.ny; j++)
   for(int i=0; i<dim.nx; i++)
   {
      v = k*dim.np + j*dim.nx + i;

      if(msk[v]<=0 || msk[v]==100 || im1[v]<=0 || im2[v]<=0)
      {
         flg[v]=0;
         continue;
      }

      im1sph.set(sim1,i,j,k);
      standardize(im1sph.v, im1sph.v, im1sph.n);

      PT[3*c + 0] = QT[3*c + 0] = i;
      PT[3*c + 1] = QT[3*c + 1] = j;
      PT[3*c + 2] = QT[3*c + 2] = k;
      detect_lm(searchsph, im2sph, sim2, PT+3*c, im1sph, QT+3*c);

      flg[v]=1;
      c++;
   }

   int n=c;

   P = (float *)calloc(3*n, sizeof(float));
   Q = (float *)calloc(3*n, sizeof(float));

   c=0;
   for(int i=0; i<dim.nv; i++)
   if(flg[i])
   {
      P[0*n + c]=(PT[3*c + 0] - (dim.nx-1)/2.0)*dim.dx;
      P[1*n + c]=(PT[3*c + 1] - (dim.ny-1)/2.0)*dim.dy;
      P[2*n + c]=(PT[3*c + 2] - (dim.nz-1)/2.0)*dim.dz;

      Q[0*n + c]=(QT[3*c + 0] - (dim.nx-1)/2.0)*dim.dx;
      Q[1*n + c]=(QT[3*c + 1] - (dim.ny-1)/2.0)*dim.dy;
      Q[2*n + c]=(QT[3*c + 2] - (dim.nz-1)/2.0)*dim.dz;

      c++;
   }

   leastSquaresAffineTrans(P, Q, n, A);

   free(PT); 
   free(QT);
   free(flg);
   free(P); 
   free(Q);
}

void generateMultiResolution(short *im, DIM dim, float *T, float *Xwarp, float *Ywarp, float *Zwarp, short *im1, short *im2, short *im4, short *im8)
{
   short *tmp;

   float I[16];
   float *Xw,*Yw,*Zw;

   Xw = (float *)calloc(dim1.nv, sizeof(float));
   Yw = (float *)calloc(dim1.nv, sizeof(float));
   Zw = (float *)calloc(dim1.nv, sizeof(float));

   for(int i=0; i<dim1.nv; i++) {
      Xw[i] = Xwarp[i];
      Yw[i] = Ywarp[i];
      Zw[i] = Zwarp[i];
   }

   combine_warps_and_trans(dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim1.dy, dim1.dz, Xw, Yw, Zw, T);

   tmp=computeReslicedImage2(im, dim, dim1, Xw, Yw, Zw);
   for(int i=0; i<dim1.nv; i++) im1[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im1, dim1, dim2, I, LIN); 
   for(int i=0; i<dim2.nv; i++) im2[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im2, dim2, dim4, I, LIN); 
   for(int i=0; i<dim4.nv; i++) im4[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im4, dim4, dim8, I, LIN); 
   for(int i=0; i<dim8.nv; i++) im8[i]=tmp[i];
   free(tmp);

   free(Xw); free(Yw); free(Zw);
}

void generateMultiResolution(short *im, DIM im_dim, float *T, short *im1, short *im2, short *im4, short *im8)
{
   short *tmp;
   float I[16];
   float *invT;		

   invT=inv4(T); 
   tmp = resliceImage(im, im_dim, dim1, invT, LIN); 
   free(invT);
   for(int i=0; i<dim1.nv; i++) im1[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im1, dim1, dim2, I, LIN); 
   for(int i=0; i<dim2.nv; i++) im2[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im2, dim2, dim4, I, LIN); 
   for(int i=0; i<dim4.nv; i++) im4[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im4, dim4, dim8, I, LIN); 
   for(int i=0; i<dim8.nv; i++) im8[i]=tmp[i];
   free(tmp);

   return;
}

void combine_warps_and_trans(float *T, float *Xwarp, float *Ywarp, float *Zwarp, float *Xout, float *Yout, float *Zout, DIM dim)
{
   int v;
   float xc,yc,zc;
   float xc1,yc1,zc1;
   float x,y,z;   
   float x1,y1,z1;   
   float i1,j1,k1;   

   xc=dim.dx*(dim.nx-1.0)/2.0;     /* +---+---+ */
   yc=dim.dy*(dim.ny-1.0)/2.0;
   zc=dim.dz*(dim.nz-1.0)/2.0;

   xc1=dim1.dx*(dim1.nx-1.0)/2.0;     /* +---+---+ */
   yc1=dim1.dy*(dim1.ny-1.0)/2.0;
   zc1=dim1.dz*(dim1.nz-1.0)/2.0;
 
   for(int k=0;k<dim.nz;k++) 
   for(int j=0;j<dim.ny;j++) 
   for(int i=0;i<dim.nx;i++) 
   {
      v = k*dim.np + j*dim.nx + i;

      // (i*dx-xc) converts from image coordinates (i,j,z) to (x,y,z) coordinates
      x = (i*dim.dx - xc);
      y = (j*dim.dy - yc);
      z = (k*dim.dz - zc);

      x1 =  T[0]*x +T[1]*y +T[2]*z  +T[3];
      y1 =  T[4]*x +T[5]*y +T[6]*z  +T[7];
      z1 =  T[8]*x +T[9]*y +T[10]*z +T[11];

      i1 = (x1 + xc1) / dim1.dx;
      j1 = (y1 + yc1) / dim1.dy;
      k1 = (z1 + zc1) / dim1.dz;

      x1 += linearInterpolator(i1,j1,k1,Xwarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      y1 += linearInterpolator(i1,j1,k1,Ywarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      z1 += linearInterpolator(i1,j1,k1,Zwarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      
      Xout[v] = x1 - i*dim.dx + xc;
      Yout[v] = y1 - j*dim.dy + yc;
      Zout[v] = z1 - k*dim.dz + zc;
   }
}

void approximate_affine(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
	int np;
	int q,N;
	float xc, yc, zc;
	float rx,ry,rz;
	float sx,sy,sz;
  	float  x,y,z;   
	float *AFF, *invAFF;	// affine transform

	double Mrx, Mry, Mrz;
	double Msx, Msy, Msz;
	double SR[9], RR[9];
	double *invRR;
	double A[9],B[3];

   AFF = (float *)calloc(16,sizeof(float));
   if(AFF==NULL)
   {
      printf("\nMemory allocation error (AFF) ...\n");
      exit(1);
   }

	np = nx*ny;

   	xc=dx*(nx-1)/2.0;     /* +---+---+ */
	yc=dy*(ny-1)/2.0;
	zc=dz*(nz-1)/2.0;

	Mrx=Mry=Mrz=0.0;
	Msx=Msy=Msz=0.0;

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes averages of the s and r vectors defined in
	// the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	N = 0;
	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		if(Xwarp[q]!=0.0 && Ywarp[q]!=0.0 && Zwarp[q]!=0.0)
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = x;
			sy = y;
			sz = z;

			Mrx += rx; Mry += ry; Mrz += rz;
			Msx += sx; Msy += sy; Msz += sz;

			N++;
		}
	}

	if(N!=0)
	{
		Mrx /= N; Mry /= N; Mrz /= N;
		Msx /= N; Msy /= N; Msz /= N;
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes the two 3x3 matrix in Eq. (2) of the
	// tech. notes.
	/////////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<9; i++) SR[i]=RR[i]=0.0;

	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		if(Xwarp[q]!=0.0 && Ywarp[q]!=0.0 && Zwarp[q]!=0.0)
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = x;
			sy = y;
			sz = z;

			rx -= Mrx; ry -= Mry; rz -= Mrz;
			sx -= Msx; sy -= Msy; sz -= Msz;

			SR[0]+=sx*rx; SR[1]+=sx*ry; SR[2]+=sx*rz;
			SR[3]+=sy*rx; SR[4]+=sy*ry; SR[5]+=sy*rz;
			SR[6]+=sz*rx; SR[7]+=sz*ry; SR[8]+=sz*rz;

			RR[0]+=rx*rx; RR[1]+=rx*ry; RR[2]+=rx*rz;
			RR[3]+=ry*rx; RR[4]+=ry*ry; RR[5]+=ry*rz;
			RR[6]+=rz*rx; RR[7]+=rz*ry; RR[8]+=rz*rz;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// estimate A according to Eq. (2) of the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	invRR = inv3(RR);
	multi(SR,3,3,invRR,3,3,A);
	free(invRR);
	/////////////////////////////////////////////////////////////////////////////////

	// estimate B according to Eq. (1) of the technical notes
	B[0] = Msx - A[0]*Mrx - A[1]*Mry - A[2]*Mrz;
	B[1] = Msy - A[3]*Mrx - A[4]*Mry - A[5]*Mrz;
	B[2] = Msz - A[6]*Mrx - A[7]*Mry - A[8]*Mrz;

	// Eq. (3) of tech. notes
	AFF[0]=(float)A[0]; AFF[1]=(float)A[1];  AFF[2]=(float)A[2];  AFF[3]=(float)B[0];
	AFF[4]=(float)A[3]; AFF[5]=(float)A[4];  AFF[6]=(float)A[5];  AFF[7]=(float)B[1];
	AFF[8]=(float)A[6]; AFF[9]=(float)A[7]; AFF[10]=(float)A[8]; AFF[11]=(float)B[2];
	AFF[12]=0.0; AFF[13]=0.0; AFF[14]=0.0; AFF[15]=1.0;

	// Eq. (3.5) of tech. notes
	invAFF = inv4(AFF);
	delete AFF;

	///////////////////////////////////////////////////////////////////////////////

	// replace T with invAFF
	for(int i=0; i<16; i++) T[i] = invAFF[i];

	delete invAFF;
}

void fillzerovectors(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
  	float  x,y,z;   
  	float  xx,yy,zz;   
	float xc,yc,zc;
	float *invT;		
	int q, np;

    np = nx*ny;

	invT=inv4(T);

   	xc=dx*(nx-1)/2.0;     /* +---+---+ */
	yc=dy*(ny-1)/2.0;
	zc=dz*(nz-1)/2.0;

	q=0;
	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		if(Xwarp[q]==0.0 && Ywarp[q]==0.0 && Zwarp[q]==0.0)
		{
		   xx = (i*dx - xc);
		   yy = (j*dy - yc);
		   zz = (k*dz - zc);

		   x = ( invT[0]*xx +invT[1]*yy +invT[2]*zz  +invT[3]  );
		   y = ( invT[4]*xx +invT[5]*yy +invT[6]*zz  +invT[7]  );
		   z = ( invT[8]*xx +invT[9]*yy +invT[10]*zz +invT[11] );

		   Xwarp[q] = x - xx;
		   Ywarp[q] = y - yy;
		   Zwarp[q] = z - zz;
	   }
	}

	free(invT);
}

void ivf(float *Xwarp1, float *Ywarp1, float *Zwarp1, DIM dim1, float *Xwarp2, float *Ywarp2, float *Zwarp2, DIM dim2)
{
   int i2, j2, k2;
   int ii, jj, kk;
   int np1, nv1;
   int np2, nv2;
   int wx, wy, wz;
   int v1, v1_part1, v1_part2;
   int v2, v2_part1, v2_part2;

   float sigma, K;
   float *W, w;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float x,y,z;
   float xx, yy, zz;
   float x2, y2, z2;
   float x1, y1, z1;

   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;

   ////////////////////////////////////////////////////////////////////////////////

   nx1 = dim1.nx; ny1 = dim1.ny; nz1 = dim1.nz;
   nx2 = dim2.nx; ny2 = dim2.ny; nz2 = dim2.nz;

   dx1 = dim1.dx; dy1 = dim1.dy; dz1 = dim1.dz;
   dx2 = dim2.dx; dy2 = dim2.dy; dz2 = dim2.dz;

   ////////////////////////////////////////////////////////////////////////////////

   nv1 = nx1 * ny1 * nz1;
   np1 = nx1 * ny1;

   nv2 = nx2 * ny2 * nz2;
   np2 = nx2 * ny2;

   ////////////////////////////////////////////////////////////////////////////////

   sigma = FWHM/2.35482; // converts FWHM to standard deviation
   K = 2.0*sigma*sigma;

   // Twice the standard deviation is greater than 95% of the area
   wx = (int)ceilf(2.0*sigma/dx2);
   wy = (int)ceilf(2.0*sigma/dy2);
   wz = (int)ceilf(2.0*sigma/dz2);

   ////////////////////////////////////////////////////////////////////////////////

//   if(opt_v)
//   {
//      printf("Standard deviation of the Gaussian kernel = %f\n\n",sigma);
//      printf("Support region size: %d x %d x %d voxels\n\n", 2*wx+1, 2*wy+1, 2*wz+1);

//      printf("Input matrix size = %d x %d x %d\n",nx1, ny1, nz1);
//      printf("Input voxel size = %7.5f x %7.5f x %7.5f mm^3\n\n",dx1, dy1, dz1);

//      printf("Output matrix size = %d x %d x %d\n",nx2, ny2, nz2);
//      printf("Output voxel size = %7.5f x %7.5f x %7.5f mm^3\n",dx2, dy2, dz2);
//   }

   ////////////////////////////////////////////////////////////////////////////////

   W = (float *)calloc(nv2,sizeof(float));
   if(W==NULL)
   {
      printf("\nMemory allocation error (W) ...\n");
      exit(1);
   }

   xc1 = dx1 * (nx1-1.0)/2.0;
   yc1 = dy1 * (ny1-1.0)/2.0;
   zc1 = dz1 * (nz1-1.0)/2.0;

   xc2 = dx2 * (nx2-1.0)/2.0;
   yc2 = dy2 * (ny2-1.0)/2.0;
   zc2 = dz2 * (nz2-1.0)/2.0;

   for(int k1=0; k1<nz1; k1++)
   {
      v1_part1 = k1*np1;
      z1 = dz1 * k1 - zc1;
      
      for(int j1=0; j1<ny1; j1++)
      {
         v1_part2 = j1*nx1;
         y1 = dy1 * j1 - yc1;

         for(int i1=0; i1<nx1; i1++)
         {
            v1 = v1_part1 + v1_part2 + i1;
            x1 = dx1 * i1 - xc1;

            if( Xwarp1[v1]!=0.0 && Ywarp1[v1]!=0.0 && Zwarp1[v1]!=0.0 )
            {
               x2 = x1 + Xwarp1[v1];
               y2 = y1 + Ywarp1[v1];
               z2 = z1 + Zwarp1[v1];

               i2 = (int)nearbyintf( (x2 + xc2)/dx2 );
               j2 = (int)nearbyintf( (y2 + yc2)/dy2 );
               k2 = (int)nearbyintf( (z2 + zc2)/dz2 );

               for(int k=-wz; k<=wz; k++)
               {
                  kk = k+k2;
                  v2_part1 = kk*np2;

                  for(int j=-wy; j<=wy; j++)
                  {
                     jj = j+j2;
                     v2_part2 = jj*nx2;

                     for(int i=-wx; i<=wx; i++)
                     {
                        ii = i+i2;

                        if( ii>0 && jj>0 && kk>0 && ii<nx2 && jj<ny2 && kk<nz2 )
                        {
                           x = dx2 * ii - xc2;
                           y = dy2 * jj - yc2;
                           z = dz2 * kk - zc2;

                           v2 = v2_part1 + v2_part2 + ii;

                           xx = x-x2;
                           yy = y-y2;
                           zz = z-z2;

                           w = expf( -(xx*xx + yy*yy + zz*zz)/K);

                           W[v2] += w;

                           Xwarp2[v2] += w*( x1 - x);
                           Ywarp2[v2] += w*( y1 - y);
                           Zwarp2[v2] += w*( z1 - z);
                        }
                     }
                  }
               }
            }
         }
      }
   }

   for(int i=0; i<nv2; i++)
   if(W[i]>0.0)
   {
      Xwarp2[i] /= W[i];
      Ywarp2[i] /= W[i];
      Zwarp2[i] /= W[i];
   }

   {
      float T[16];
      approximate_affine(nx2, ny2, nz2, dx2, dy2, dz2, Xwarp2, Ywarp2, Zwarp2, T);
      fillzerovectors(nx2, ny2, nz2, dx2, dy2, dz2, Xwarp2, Ywarp2, Zwarp2, T);
   }

   free(W);
}

void read_default_atlas_names(const char *brainwashatlasdir, int &natlas)
{
   int2 L; // number of charaters in the filename (i.e., its length)

   char fileprefix[1024];
   char filename[1024];

   DIR *dp;
   // struct direct *dir; // this one didn't work on SUN 
   struct dirent *dir; 

   natlas=0;

   // open the directory and associate the directory stream dp with it 
   dp=opendir(brainwashatlasdir);

   // If opendir() returns a NULL pointer, then brainwashatlasdir cannot be accessed, 
   // or if it cannot malloc() enough memory to hold the whole thing. 
   if(dp==NULL)
   {
      printf("Warning: Cannot access %s, returning ...\n",brainwashatlasdir);
      return;
   }

   // Get a pointer to the next directory entry using readdir().  If readdir()
   // returns NULL, the end of the directory has been reached.  
   while( (dir=readdir(dp)) != NULL && natlas<MAXNATLAS)
   {
      L=strlen(dir->d_name);

      if(L != 8 || strcmp(dir->d_name+(L-4), ".nii")!=0 )
         continue;

      sprintf(filename,"%s/%s",brainwashatlasdir,dir->d_name);

      // then looks like a real NIFTI image since niftiFilename does some checks
      if( niftiFilename(fileprefix, filename)==1 ) 
      {
         sprintf(atlasfilename[natlas],"%s",fileprefix);
         natlas++; // increment natlas to indicate that one more atlas has been read 
      }
   };

   // close the directory stream and free the structure associated with dp 
   closedir(dp);
}

void read_atlas_list(const char *atlaslistfile, int &natlas)
{
   FILE *fp;
   char filename[1024];
   char fileprefix[1024];

   natlas=0;

   fp=fopen(atlaslistfile,"r");

   if(fp==NULL)
   {
      printf("Warning: Cannot access %s, returning ...\n",atlaslistfile);
      return;
   }

   while( fscanf(fp,"%s",filename) != EOF && natlas<MAXNATLAS)
   {
      if(strlen(filename) < 5 )
         continue;

      // then looks like a real NIFTI image since niftiFilename does some checks
      if( niftiFilename(fileprefix, filename)==1 ) 
      {
         sprintf(atlasfilename[natlas],"%s",filename);
         natlas++; // increment natlas to indicate that one more atlas has been read 
      }
   };

   fclose(fp);
}

void matchpatch(SPH &searchsph, SPH &testsph, SHORTIM testim, int P[], SPH &refsph)
{
   for(int n=0; n<searchsph.n; n++)
   {
      testsph.set(testim, P[0]+searchsph.i[n], P[1]+searchsph.j[n], P[2]+searchsph.k[n]);
      standardize(testsph.v, testsph.n);

      searchsph.v[n] = dot(testsph.v, refsph.v, testsph.n);
   }
}

int main(int argc, char **argv)
{
   char orient[4]="";
   float *invT;		
   short *tmpmsk;
   char atlaspath[1024];
   char atlasmskpath[1024];
   int natlas_used=11;

   int4 *atlas_indx;

   int2 *PILbraincloud;
   DIM PILbraincloudDim;
   nifti_1_header PILbraincloudHdr;

   float threshold=50.0;

   char atlaslistfile[1024]=""; 

   float *label, *evidence0, *evidence1;

   char brainwashatlasdir[1024]=""; 
   int natlas; // number of available atlases
   
   // input image full path
   char subImageFile[1024]=""; 


   // extracted image filename without suffix
   char subprefix[1024]=""; 
   char atlprefix[1024]="";

   // original atlas and subject images loaded 
   short *sub, *atl;

   // NIFTI headers for subject and atlas images
   nifti_1_header subHdr, atlHdr;

   // structures for holding atlas and subject image dimensions
   DIM subDim, atlDim;

   // a linear transformation that brings the atlas/subject image from its native space to PIL space
   float sub2PIL[16];
   float *atl2PIL;

   // filenames where AC, PC, and VSPS are manually specific for atlas/subject image
   char subLMfile[1024]="";

   // number if iterations of NL registration at each resolution level
   int iter8=4;
   int iter4=2;
   int iter2=1;
   int iter1=1;

   // Subject and atlas PIL versions and brain mask at different resolutons
   short *PILsub1=NULL, *PILsub2=NULL, *PILsub4=NULL, *PILsub8=NULL;
   short *PILatl1=NULL, *PILatl2=NULL, *PILatl4=NULL, *PILatl8=NULL;
   short *PILmsk1=NULL, *PILmsk2=NULL, *PILmsk4=NULL, *PILmsk8=NULL;

   // intermediate displacement field from PILatl to PILsub
   float *Xwarp=NULL, *Ywarp=NULL, *Zwarp=NULL;

   // 4x4 identity matrix
   float I[16];

   // space holder for a generic filename
   char filename[1024];

   int patch_r=3; // patch radius
   int search_R=3; // search radius 

   // atlas to subject image affine transformation
   float *atl_to_sub;

   opt_ppm = NO;
   opt_txt=NO;

   /////////////////////////////////////////////////////////////////////
   // maximize the stack size
   /////////////////////////////////////////////////////////////////////
   {
      // define a rlimit structure
      struct rlimit rlpty;

      // get the resource limit of type RLIMIT_STACK for stack size
      (void) getrlimit(RLIMIT_STACK,&rlpty);

      // set the stack size to its maximum
      rlpty.rlim_cur=rlpty.rlim_max;
      (void) setrlimit(RLIMIT_STACK,&rlpty);
   }
   /////////////////////////////////////////////////////////////////////

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'V':
            printf("Version 2.0 (last modified: May 17, 2016)\n");
            exit(0);
         case 't':
            threshold = atof(optarg);
            break;
         case 'r':
            patch_r = atoi(optarg);
            break;
         case 'p':
            opt_ppm=YES;
            break;
         case 'R':
            search_R = atoi(optarg);
            break;
         case 'n':
            natlas_used = atoi(optarg);
            break;
         case '8':
            iter8 = atoi(optarg);
            break;
         case '4':
            iter4 = atoi(optarg);
            break;
         case '2':
            iter2 = atoi(optarg);
            break;
         case '1':
            iter1 = atoi(optarg);
            break;
         case 'l':
            sprintf(atlaslistfile,"%s",optarg);
            break;
         case 'd':
            sprintf(brainwashatlasdir,"%s",optarg);
            break;
         case 'M':
            sprintf(subLMfile,"%s",optarg);
            break;
         case 'h':
            print_help_and_exit();
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'i':
            sprintf(subImageFile,"%s",optarg);
            break;
         case '?':
            print_help_and_exit();
		}
	}

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Variable initializations

   getARTHOME();

   // revert to the default value if unreasonable input
   if(threshold<0.0 || threshold>100.0) threshold=50.0;

   // Initialize dim1, dim2, dim4, and dim8 multi-resolution image dimensions
   dim1.nx = XMATRIXSIZE;
   dim1.ny = YMATRIXSIZE;
   dim1.nz = ZMATRIXSIZE;
   dim1.np=dim1.nx*dim1.ny; 
   dim1.nv=dim1.np*dim1.nz; 
   dim1.dx = dim1.dy = dim1.dz = VOXELSIZE;
   PILatl1 = (short *)calloc(dim1.nv,sizeof(short));
   PILsub1 = (short *)calloc(dim1.nv,sizeof(short));
   PILmsk1 = (short *)calloc(dim1.nv,sizeof(short));

   dim2.nx = XMATRIXSIZE/2;
   dim2.ny = YMATRIXSIZE/2;
   dim2.nz = ZMATRIXSIZE/2;
   dim2.np=dim2.nx*dim2.ny; 
   dim2.nv=dim2.np*dim2.nz; 
   dim2.dx = dim2.dy = dim2.dz = VOXELSIZE*2.0;
   PILatl2 = (short *)calloc(dim2.nv,sizeof(short));
   PILsub2 = (short *)calloc(dim2.nv,sizeof(short));
   PILmsk2 = (short *)calloc(dim2.nv,sizeof(short));

   dim4.nx = XMATRIXSIZE/4;
   dim4.ny = YMATRIXSIZE/4;
   dim4.nz = ZMATRIXSIZE/4;
   dim4.np=dim4.nx*dim4.ny; 
   dim4.nv=dim4.np*dim4.nz; 
   dim4.dx = dim4.dy = dim4.dz = VOXELSIZE*4.0;
   PILatl4 = (short *)calloc(dim4.nv,sizeof(short));
   PILsub4 = (short *)calloc(dim4.nv,sizeof(short));
   PILmsk4 = (short *)calloc(dim4.nv,sizeof(short));

   dim8.nx = XMATRIXSIZE/8;
   dim8.ny = YMATRIXSIZE/8;
   dim8.nz = ZMATRIXSIZE/8;
   dim8.np=dim8.nx*dim8.nx; 
   dim8.nv=dim8.np*dim8.nz; 
   dim8.dx = dim8.dy = dim8.dz = VOXELSIZE*8.0;
   PILatl8 = (short *)calloc(dim8.nv,sizeof(short));
   PILsub8 = (short *)calloc(dim8.nv,sizeof(short));
   PILmsk8 = (short *)calloc(dim8.nv,sizeof(short));

   // reset unreasonable values to default
   if(iter8 < 0 || iter8>MAXITER ) iter8=4;
   if(iter4 < 0 || iter4>MAXITER ) iter4=3;
   if(iter2 < 0 || iter2>MAXITER ) iter2=2;
   if(iter1 < 0 || iter1>MAXITER ) iter1=0;

   if(patch_r <= 0 || patch_r > MAXR ) patch_r=3;
   if(search_R <= 0 || search_R > MAXR ) search_R=3;

   if(brainwashatlasdir[0]=='\0')
   {
      sprintf(brainwashatlasdir,"%s/brainwashatlas",ARTHOME);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////

   if(atlaslistfile[0]=='\0')
   {
      if(opt_v) printf("Brainwash atlases are read from %s directory.\n",brainwashatlasdir);
      read_default_atlas_names(brainwashatlasdir, natlas);
   }
//   else
//   {
//      if(opt_v) printf("Brainwash atlases are read from %s list.\n",atlaslistfile);
//      read_atlas_list(atlaslistfile, natlas);
//   }

   if(opt_v) 
   {
      printf("Number of available atlases = %d \n",natlas);
   }

   if(natlas<=0)
   {
      exit(0);
   }

   atlas_indx = (int *)calloc(natlas, sizeof(int));

   ////////////////////////////////////////////////////////////////////////////////////////////

   // compute PILmsk1, PILmsk2, PILmsk4 and PILmsk8
   {
      sprintf(filename,"%s/PILbrain.nii",ARTHOME);

      PILbraincloud = (int2 *)read_nifti_image(filename, &PILbraincloudHdr);

      if(PILbraincloud==NULL)
      {
            printf("Error reading %s, aborting ...\n", filename);
            exit(1);
      }

      set_dim(PILbraincloudDim, PILbraincloudHdr);

      set_to_I(I,4);
      generateMultiResolution(PILbraincloud, PILbraincloudDim, I, PILmsk1, PILmsk2, PILmsk4, PILmsk8);
   }

   // ensure the subject image is specified at the command line
   if( subImageFile[0] == '\0')
   {
      printf("\nPlease specify an image for skull-stripping using -i <image.nii> ...\n\n");		
      exit(0);
   } 

   // extract the subject filename without path/suffix
   if( niftiFilename(subprefix, subImageFile)==0 ) { exit(0); }

   // read the subject image
   sub = (int2 *)read_nifti_image(subImageFile, &subHdr);
   if(sub==NULL) 
   {
      printf("\nError: Reading image %s failed.\n\n",subImageFile);
      exit(0);
   }
   set_dim(subDim,subHdr); // transfer info from subHdr to subDim

   // print some info about the image
   if(opt_v)
   {
      printf("Image file: %s\n",subImageFile);
      printf("\tMatrix size: %d x %d x %d\n",subDim.nx, subDim.ny, subDim.nz);
      printf("\tVoxel size: %6.4f x %6.4f x %6.4f\n",subDim.dx, subDim.dy, subDim.dz);
   }

   // find sub2PIL using automated MSP, AC/PC and 8 MSP landmarks
   if(opt_v) printf("Computing PIL transformation for %s ...\n",subImageFile);
   opt_ppm=YES;
   if(subLMfile[0] != '\0' && opt_v) printf("Image landmarks are read from %s\n",subLMfile);
   orient[0]='\0';
   new_PIL_transform(subImageFile, subLMfile, orient, sub2PIL, 1);
   opt_ppm=NO;

   // atlas selection
   {
      if(opt_v) printf("Selecting top %d atlases ...\n",natlas_used);

      float4 *corr;

      int2 *tmp1, *tmp2;

      corr=(float *)calloc(natlas,sizeof(float));
      tmpmsk =(int2 *)calloc(PILbraincloudDim.nv,sizeof(int2));

      for(int v=0; v<PILbraincloudDim.nv; v++) 
      {
         if(PILbraincloud[v]<100 && PILbraincloud[v]>0) 
         {
            tmpmsk[v]=1; 
         } else {
            tmpmsk[v]=0;
         }
      }
      //save_nifti_image("M1.nii", tmpmsk, &PILbraincloudHdr);

      invT=inv4(sub2PIL); 
      tmp1 = resliceImage(sub, subDim, PILbraincloudDim, invT, LIN); 
      free(invT);

      for(int a=0; a<natlas; a++)
      {
         sprintf(atlaspath,"%s/%s_PIL.nii",brainwashatlasdir,atlasfilename[a]);

         tmp2 = (int2 *)read_nifti_image(atlaspath, &PILbraincloudHdr);

         atlas_indx[a]=a;
         corr[a] = pearsonCorrelation(tmp1, tmp2, tmpmsk, PILbraincloudDim.nv);

         free(tmp2);
      }

      hpsort(natlas, corr, atlas_indx);

      if(opt_v)
      {
         printf("Selected atlases:\n");
         for(int i=0; i<natlas_used; i++)
         {
            int a = atlas_indx[natlas-1-i];
            float c=corr[natlas-1-i];
            printf("Rank=%03d, Atlas=%s, Correlation=%5.3f\n",i+1,atlasfilename[a],c);
         }
      }

      free(corr);
      free(tmp1);
      free(tmpmsk);
   }

   if(opt_v) printf("Patch radius = %d mm\n", patch_r);
   if(opt_v) printf("Search radius = %d mm\n", search_R);

   label = (float *)calloc(dim1.nv,sizeof(float));
   evidence0 = (float *)calloc(dim1.nv,sizeof(float));
   evidence1 = (float *)calloc(dim1.nv,sizeof(float));

   // compute PILsub1, PILsub2, PILsub4 and PILsub8
   generateMultiResolution(sub, subDim, sub2PIL, PILsub1, PILsub2, PILsub4, PILsub8);

   atl_to_sub=(float *)calloc(natlas_used*16,sizeof(float));
   atl2PIL=(float *)calloc(natlas_used*16,sizeof(float));

   for(int i=0; i<natlas_used; i++)
   {
      int a; // atlas index
      short *atlmsk;
      float A[16];
      short *tmp;

      a = atlas_indx[natlas-2-i];

      sprintf(atlaspath,"%s/%s.nii",brainwashatlasdir,atlasfilename[a]);
      sprintf(atlasmskpath,"%s/%s_msk.nii",brainwashatlasdir,atlasfilename[a]);

      // extract the atlas filename without path/suffix
      if( niftiFilename(atlprefix, atlasmskpath )==0 ) { exit(0); } // just for checking
      if( niftiFilename(atlprefix, atlaspath )==0 ) { exit(0); }

      // read the atlas image
      atl = (int2 *)read_nifti_image(atlaspath, &atlHdr);
      if(atl==NULL) 
      {
         printf("\nError: Reading atlas %s failed.\n\n",atlaspath);
         exit(0);
      }
      set_dim(atlDim,atlHdr); // transfer info from atlHdr to atlDim

      atlmsk = (int2 *)read_nifti_image(atlasmskpath, &atlHdr);
      if(atlmsk==NULL) 
      {
         printf("\nError: Reading atlas %s failed.\n\n",atlasmskpath);
         exit(0);
      }

      for(int v=0; v<atlDim.nv; v++) if(atlmsk[v]==0) atl[v]=0;

      // print some info about the atlas image
      if(opt_v)
      {
         printf("Atlas %03d: %s\n",i+1,atlaspath);
      }

      if(opt_v) printf("Computing PIL transformation for %s ...\n", atlaspath);
      orient[0]='\0';
      new_PIL_transform(atlaspath,"", orient, atl2PIL+i*16, 1);

      ////////////////////////////////////////////////////////////////////////////////////////////
      // Compute the affine transformation  atl_to_sub and update atl2PIL
      if(opt_v) printf("Affine registration @ 12.5%% resolution ...\n");
      generateMultiResolution(atl, atlDim, atl2PIL+i*16, PILatl1, PILatl2, PILatl4, PILatl8);
      affineReg(PILatl8, PILsub8, PILmsk8, dim8, patch_r, search_R, A);
      multi(A,4,4,atl2PIL+i*16,4,4,atl2PIL+i*16);  // update atl2PIL

      if(opt_v) printf("Affine registration @ 25%% resolution ...\n");
      generateMultiResolution(atl, atlDim, atl2PIL+i*16, PILatl1, PILatl2, PILatl4, PILatl8);
      affineReg(PILatl4, PILsub4, PILmsk4, dim4, patch_r, search_R, A);
      multi(A,4,4,atl2PIL+i*16,4,4,atl2PIL+i*16);  // update atl2PIL

      if(opt_v) printf("Affine registration @ 50%% resolution ...\n");
      generateMultiResolution(atl, atlDim, atl2PIL+i*16, PILatl1, PILatl2, PILatl4, PILatl8);
      affineReg(PILatl2, PILsub2, PILmsk2, dim2, patch_r, search_R, A);
      multi(A,4,4,atl2PIL+i*16,4,4,atl2PIL+i*16);  // update atl2PIL

      invT = inv4(sub2PIL);
      multi(invT,4,4, atl2PIL+i*16,4,4, atl_to_sub+i*16);
      free(invT);

      if(opt_v) printMatrix(atl_to_sub+i*16,4,4,"Atlas image -> subject image affine transformation",NULL);

      invT=inv4(atl2PIL+i*16); 
      tmp = resliceImage(atlmsk, atlDim, dim1, invT, LIN); 
      free(invT);
      for(int v=0; v<dim1.nv; v++) label[v] += tmp[v]/100.0;
      free(tmp);

      free(atlmsk);
      free(atl);
   }

   tmpmsk = (short *)calloc(dim1.nv,sizeof(short));
   for(int v=0; v<dim1.nv; v++)
   {
      label[v] /= natlas_used; // makes label range from 0.0-1.0

      if(label[v]>0.0 && label[v]<1.0) { tmpmsk[v]=1; } else { tmpmsk[v]=0; }

   }
   //save_nifti_image("M2.nii", tmpmsk, &PILbraincloudHdr);
   
   for(int i=0; i<dim1.nv; i++) evidence0[i]=evidence1[i]=0.0; //resetting

   for(int ai=0; ai<natlas_used; ai++)
   {
      int a; // atlas index
      short *atlmsk;
      int v, vv; // voxel index
      SPH subsph(patch_r);
      SPH atlsph(patch_r);
      SPH searchsph(search_R);
      int P[3];
      int Q[3];
      float cc;

      a = atlas_indx[natlas-2-ai];

      sprintf(atlaspath,"%s/%s.nii",brainwashatlasdir,atlasfilename[a]);
      sprintf(atlasmskpath,"%s/%s_msk.nii",brainwashatlasdir,atlasfilename[a]);

      // extract the atlas filename without path/suffix
      if( niftiFilename(atlprefix, atlasmskpath )==0 ) { exit(0); } // just for checking
      if( niftiFilename(atlprefix, atlaspath )==0 ) { exit(0); }

      // read the atlas image
      atl = (int2 *)read_nifti_image(atlaspath, &atlHdr);
      if(atl==NULL) 
      {
         printf("\nError: Reading atlas %s failed.\n\n",atlaspath);
         exit(0);
      }
      set_dim(atlDim,atlHdr); // transfer info from atlHdr to atlDim

      atlmsk = (int2 *)read_nifti_image(atlasmskpath, &atlHdr);
      if(atlmsk==NULL) 
      {
         printf("\nError: Reading atlas %s failed.\n\n",atlasmskpath);
         exit(0);
      }

      // print some info about the atlas image
      if(opt_v)
      {
         printf("Atlas %03d: %s\n",ai+1,atlaspath);
      }

      if(PILatl1!=NULL) free(PILatl1);
      if(PILmsk1!=NULL) free(PILmsk1);

      invT=inv4(atl2PIL+ai*16); 
      PILmsk1 = resliceImage(atlmsk, atlDim, dim1, invT, NEARN); 
      free(invT);

      invT=inv4(atl2PIL+ai*16); 
      PILatl1 = resliceImage(atl, atlDim, dim1, invT, LIN); 
      free(invT);

      //sprintf(filename,"%s_bw.nii",atlprefix);
      //save_nifti_image(filename, PILatl1, &PILbraincloudHdr);

      SHORTIM subim; 
      set_dim(subim,dim1); 
      subim.v=PILsub1;

      SHORTIM atlim; 
      set_dim(atlim,dim1); 
      atlim.v=PILatl1;

      for(int k=0; k<dim1.nz; k++)
      for(int j=0; j<dim1.ny; j++)
      for(int i=0; i<dim1.nx; i++)
      {
         v = k*dim1.np + j*dim1.nx + i;
         if( tmpmsk[v] == 0 ) continue;
         
         subsph.set(subim,i,j,k);
         standardize(subsph.v,subsph.n);

         P[0]=i; P[1]=j; P[2]=k;
         Q[0]=i; Q[1]=j; Q[2]=k;
         cc=0.0;
         cc=detect_lm(searchsph, atlsph, atlim, P, subsph, Q);

         vv = Q[2]*dim1.np + Q[1]*dim1.nx + Q[0];

         if( PILmsk1[vv] == 100 ) evidence1[v] += cc;
         if( PILmsk1[vv] == 0 )   evidence0[v] += cc;
      }
  
      free(PILmsk1); PILmsk1=NULL;
      free(PILatl1); PILatl1=NULL;
      free(atl);
      free(atlmsk);
   }

   for(int v=0; v<dim1.nv; v++)
   if( tmpmsk[v] != 0 )
   {
//      if( evidence1[v]*label[v] > evidence0[v]*(1.0-label[v]) ) label[v]=1.0; else label[v]=0.0;
      if( evidence1[v] > evidence0[v] ) label[v]=1.0; else label[v]=0.0;
   }

   free(tmpmsk);

   {
      short *tmp;

      PILmsk1 = (int2 *)calloc(dim1.nv, sizeof(int2));

      for(int v=0; v<dim1.nv; v++) if(label[v]>0.5) PILmsk1[v] = 100; else PILmsk1[v]=0;

      tmp = resliceImage(PILmsk1, dim1, subDim, sub2PIL, LIN); 

      for(int v=0; v<subDim.nv; v++) if(tmp[v]<50) sub[v] = 0;

      sprintf(filename,"%s_bw.nii",subprefix);
      save_nifti_image(filename, sub, &subHdr);

      free(tmp);
   }
exit(0);

   Xwarp = (float *)calloc(dim1.nv, sizeof(float));
   Ywarp = (float *)calloc(dim1.nv, sizeof(float));
   Zwarp = (float *)calloc(dim1.nv, sizeof(float));

   for(int ai=0; ai<natlas_used; ai++)
   {
      int a; // atlas index
      short *atlmsk;
      short *tmp;
      float *Xout=NULL, *Yout=NULL, *Zout=NULL; 

      a = atlas_indx[natlas-2-ai];

      sprintf(atlaspath,"%s/%s.nii",brainwashatlasdir,atlasfilename[a]);
      sprintf(atlasmskpath,"%s/%s_msk.nii",brainwashatlasdir,atlasfilename[a]);

      // extract the atlas filename without path/suffix
      if( niftiFilename(atlprefix, atlasmskpath )==0 ) { exit(0); } // just for checking
      if( niftiFilename(atlprefix, atlaspath )==0 ) { exit(0); }

      // read the atlas image
      atl = (int2 *)read_nifti_image(atlaspath, &atlHdr);
      if(atl==NULL) 
      {
         printf("\nError: Reading atlas %s failed.\n\n",atlaspath);
         exit(0);
      }
      set_dim(atlDim,atlHdr); // transfer info from atlHdr to atlDim

      atlmsk = (int2 *)read_nifti_image(atlasmskpath, &atlHdr);
      if(atlmsk==NULL) 
      {
         printf("\nError: Reading atlas %s failed.\n\n",atlasmskpath);
         exit(0);
      }

      // print some info about the atlas image
      if(opt_v)
      {
         printf("Atlas %03d: %s\n",ai+1,atlaspath);
      }

      //////////////////////////////////////////////////////////////////
      // Non-linear registration from PILatl to PILsub
      
      for(int v=0; v<dim1.nv; v++) Xwarp[v]=Ywarp[v]=Zwarp[v]=0.0;

/*
      if(opt_v) printf("Non-linear registration @ 12.5%% resolution ...\n");
      for(int i=0; i<iter8; i++)
      {
         generateMultiResolution(atl, atlDim, atl2PIL+ai*16, Xwarp, Ywarp, Zwarp, PILatl1, PILatl2, PILatl4, PILatl8);
         brainwashnlReg(PILatl8, PILsub8, PILmsk8, dim8, patch_r, search_R/3, Xwarp, Ywarp, Zwarp,5);
      }

      if(opt_v) printf("Non-linear registration @ 25%% resolution ...\n");
      for(int i=0; i<iter4; i++)
      {
         generateMultiResolution(atl, atlDim, atl2PIL+ai*16, Xwarp, Ywarp, Zwarp, PILatl1, PILatl2, PILatl4, PILatl8);
         brainwashnlReg(PILatl4, PILsub4, PILmsk4, dim4, patch_r, 2*search_R/3, Xwarp, Ywarp, Zwarp,5);
      }
*/

      if(opt_v) printf("Non-linear registration @ 50%% resolution ...\n");
      for(int i=0; i<iter2; i++)
      {
         generateMultiResolution(atl, atlDim, atl2PIL+ai*16, Xwarp, Ywarp, Zwarp, PILatl1, PILatl2, PILatl4, PILatl8);
         brainwashnlReg(PILatl2, PILsub2, PILmsk2, dim2, patch_r, search_R, Xwarp, Ywarp, Zwarp,5);
      }

      if(opt_v) printf("Non-linear registration @ 100%% resolution ...\n");
      for(int i=0; i<iter1; i++)
      {
         generateMultiResolution(atl, atlDim, atl2PIL+ai*16, Xwarp, Ywarp, Zwarp, PILatl1, PILatl2, PILatl4, PILatl8);
         brainwashnlReg(PILatl1, PILsub1, PILmsk1, dim1, patch_r, search_R, Xwarp, Ywarp, Zwarp,5);
      }

      //////////////////////////////////////////////////////////////////
      // combine atl2PIL, (Xwarp,Ywarp,Zwarp) and sub2PIL to obtain (Xout,Yout,Zout)

      Xout = (float *)calloc(subDim.nv, sizeof(float));
      Yout = (float *)calloc(subDim.nv, sizeof(float));
      Zout = (float *)calloc(subDim.nv, sizeof(float));

      combine_warps_and_trans(dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim1.dy, dim1.dz, Xwarp, Ywarp, Zwarp, atl2PIL+ai*16);
      combine_warps_and_trans(sub2PIL, Xwarp, Ywarp, Zwarp, Xout, Yout, Zout, subDim);

      tmp=computeReslicedImage2(atlmsk,atlDim,subDim,Xout,Yout,Zout);

      for(int v=0; v<subDim.nv; v++) label[v] += tmp[v]/100.0;
      free(tmp);

      free(Xout); free(Yout); free(Zout);

      free(atlmsk);
      free(atl);
   }

   for(int i=0; i<subDim.nv; i++)
   {
      label[i] /= natlas_used; // makes label range from 0.0-1.0
   }

   {
      short *tmp; // brain cloud in PIL space

      tmp = (int2 *)calloc(subDim.nv, sizeof(int2));

      for(int i=0; i<subDim.nv; i++) tmp[i] = (int2)(label[i]*100.0 + 0.5);

      for(int v=0; v<subDim.nv; v++) if( tmp[v]>threshold) tmp[v]=1; else tmp[v]=0;

      save_nifti_image("kk.nii", tmp, &subHdr);

      free(tmp);
   }

   //sprintf(filename,"%s_brainwash.nii",subprefix);
   //save_nifti_image(filename, sub, &subHdr);

/*
   for(int i=0; i<subDim.nv; i++)
      sub[i] = (short)(label[i]+0.5);
   sprintf(filename,"%s_brainmask.nii",subprefix);
   save_nifti_image(filename, sub, &subHdr);
*/

   ////////////////////////////////////////////////////////////////////////////////////////////
   // free allocated memory 
   free(PILbraincloud);
   free(Xwarp); free(Ywarp); free(Zwarp);
   free(sub);
   free(PILatl1); free(PILatl2); free(PILatl4); free(PILatl8);
   free(PILsub1); free(PILsub2); free(PILsub4); free(PILsub8);
   free(PILmsk1); free(PILmsk2); free(PILmsk4); free(PILmsk8);
   free(label);
   free(atlas_indx);
   ////////////////////////////////////////////////////////////////////////////////////////////

   if(opt_v) printf("THE END\n");
}
