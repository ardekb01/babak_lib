#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <babak_lib.h>
#include <resizeX.h>

short *resizeZ(float *image1, int nx1, int ny1, int nz1, float dz1, int nz2, float dz2);
float *resizeZ(float *image1, int nx1, int ny1, int nz1, float dz1, int *nz2, float *dz2);

float *resizeY(float *image1, int nx1, int ny1, int nz1, float dy1, int ny2, float dy2);
short *resizeY(float *image1, int nx1, int ny1, float dy1, int ny2, float dy2);
float *resizeY(float *image1, int nx1, int ny1, float dy1, int ny2, float *dy2);

float *resizeXY(float *image1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2)
{
	float *imA;
	float *image2;

	imA = resizeX(image1, nx1, ny1, dx1, nx2, dx2);

	image2 = resizeY(imA, nx2, ny1, dy1, ny2, &dy2);
	free(imA);
   
	return(image2);
}

float *resizeY(float *image1, int nx1, int ny1, float dy1, int ny2, float *dy2)
{
   int n;
   int np2;
   int j0;
   float d,dl;
   float sd=0.0;
   float *h;
   float *image2;
   float *x;
   float yc1,yc2;

   np2=nx1*ny2;

   yc1 = dy1*(ny1-1.0)/2.0;
   yc2 = (*dy2)*(ny2-1.0)/2.0;

   image2=(float *)calloc(np2,sizeof(float));
   if(image2==NULL) return(NULL);

   if(dy1 < (*dy2) )
   {
      sd=(float)sqrt( (0.5/log(2.0))*( (*dy2)*(*dy2) - dy1*dy1 ) );
      sd /= dy1;
   } else sd=0.0;

   h = gaussian_kernel(sd,&n);

   x=(float *)calloc(ny1,sizeof(float));

   for(int i=0;i<nx1;i++)
   {
      for(int l=0;l<ny1;l++)
         x[l]=image1[nx1*l +i];

      for(int j=0;j<ny2;j++)
      {
         d = (j*(*dy2) - yc2 + yc1)/dy1;
         j0=(int)d;
         dl=d - j0;

         if(dl==0.0)
            image2[nx1*j+i]=(conv_pnt_sk(x,ny1,h,n,j0));
         else
            image2[nx1*j+i]=( (1.0-dl)*conv_pnt_sk(x,ny1,h,n,j0) + dl*conv_pnt_sk(x,ny1,h,n,j0+1) );
      }
   }

   free(h);
   free(x);

   return(image2);
}

short *resizeXY(short *image1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2)
{
	float *imA;
	short *image2;

	imA = resizeX(image1, nx1, ny1, dx1, nx2, dx2);

	image2 = resizeY(imA, nx2, ny1, dy1, ny2, dy2);
	free(imA);
   
	return(image2);
}

short *resizeY(float *image1, int nx1, int ny1, float dy1, int ny2, float dy2)
{
   int n;
   int np2;
   int j0;
   float d,dl;
   float sd=0.0;
   float *h;
   short *image2;
   float *x;
   float yc1,yc2;

   np2=nx1*ny2;

   yc1 = dy1*(ny1-1.0)/2.0;
   yc2 = dy2*(ny2-1.0)/2.0;

   image2=(short *)calloc(np2,sizeof(short));
   if(image2==NULL) return(NULL);

   if(dy1 < dy2)
   {
      sd=(float)sqrt( (0.5/log(2.0))*( dy2*dy2 - dy1*dy1 ) );
      sd /= dy1;
   } else sd=0.0;

   h = gaussian_kernel(sd,&n);

   x=(float *)calloc(ny1,sizeof(float));

   for(int i=0;i<nx1;i++)
   {
      for(int l=0;l<ny1;l++)
         x[l]=image1[nx1*l +i];

      for(int j=0;j<ny2;j++)
      {
         d = (j*dy2 - yc2 + yc1)/dy1;
         j0=(int)d;
         dl=d - j0;

         if(dl==0.0)
            image2[nx1*j+i]=(short)(conv_pnt_sk(x,ny1,h,n,j0) + 0.5);
         else
            image2[nx1*j+i]=(short)( (1.0-dl)*conv_pnt_sk(x,ny1,h,n,j0) + dl*conv_pnt_sk(x,ny1,h,n,j0+1) + 0.5 );
      }
   }

   free(h);
   free(x);

   return(image2);
}

short *resizeXYZ(char *image1, 
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2)
{
	float *imA, *imB;
	short *image2;
	short *tmp;
	int nv1;

	nv1 = nx1*ny1*nz1;

	tmp = (short *)calloc(nv1,sizeof(short));
	for(int i=0; i<nv1; i++) tmp[i] = image1[i];
	imA = resizeX(tmp, nx1, ny1, nz1, dx1, nx2, dx2);
	free(tmp);

	imB = resizeY(imA, nx2, ny1, nz1, dy1, ny2, dy2);
	free(imA);
   
	image2=resizeZ(imB, nx2, ny2, nz1, dz1, nz2, dz2);
	free(imB);

	return(image2);
}

uchar *resizeZ_UC(float *image1, int nx1, int ny1, int nz1, float dz1, int nz2, float dz2)
{
	int n;
	int np1;
	int k0;
	float d,dl;
	float sd=0.0;
	float *h;
	float *x;
	float zc1,zc2;
	uchar *image2;

	np1=nx1*ny1;

	zc1 = dz1*(nz1-1.0)/2.0;
	zc2 = dz2*(nz2-1.0)/2.0;

	image2=(uchar *)calloc(np1*nz2,1);
	if(image2==NULL) return(NULL);

	if(nz2 < nz1)
	{
		sd=(float)sqrt( (0.5/log(2.0))*( dz2*dz2 - dz1*dz1 ) );

		sd /= dz1;

		h = gaussian_kernel(sd,&n);

		x=(float *)calloc(nz1,sizeof(float));

		for(int j=0;j<ny1;j++)
		for(int i=0;i<nx1;i++)
		{
			for(int l=0;l<nz1;l++)
				x[l]=image1[np1*l + nx1*j +i];

			for(int k=0;k<nz2;k++)
			{
				d = (k*dz2 - zc2 + zc1)/dz1;
				k0=(int)d;
				dl=d - k0;

				if(dl==0)
					image2[np1*k+nx1*j+i]=(uchar)(conv_pnt_sk(x,nz1,h,n,k0)+0.5);
				else
					image2[np1*k+nx1*j+i]=(uchar) ( (1.0-dl)*conv_pnt_sk(x,nz1,h,n,k0) +
                   	dl*conv_pnt_sk(x,nz1,h,n,k0+1) + 0.5 );
			}
		}

		free(h);
		free(x);
	}
	else if(nz2 > nz1)
	{
		for(int j=0;j<ny1;j++)
		for(int i=0;i<nx1;i++)
		for(int k=0;k<nz2;k++)
		{
			d = (k*dz2 - zc2 + zc1)/dz1;
			k0=(int)d;
			dl=d - k0;

			if(k0<0 || k0>=(nz1-1))
			{
				image2[np1*k+nx1*j+i]=0;
				continue;
			}

			if(dl==0.0)
				image2[np1*k+nx1*j+i]=(uchar)(image1[np1*k0 + nx1*j +i]+0.5);
			else
				image2[np1*k+nx1*j+i]=(uchar)((1.0-dl)*image1[np1*k0 + nx1*j +i]
					+ dl*image1[np1*(k0+1) + nx1*j +i] + 0.5 );
		}
	}
	else
	{
		for(int i=0;i<np1*nz1;i++)
			image2[i]=(uchar)(image1[i]+0.5);
	}

   return(image2);
}

uchar *resizeXYZ(uchar *image1, 
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2)
{
	float *imA, *imB;
	uchar *image2;

	imA = resizeX(image1, nx1, ny1, nz1, dx1, nx2, dx2);

	imB = resizeY(imA, nx2, ny1, nz1, dy1, ny2, dy2);
	free(imA);
   
	image2=resizeZ_UC(imB, nx2, ny2, nz1, dz1, nz2, dz2);
	free(imB);

	return(image2);
}

short *resizeXYZ(short *image1, 
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2)
{
	float *imA, *imB;
	short *image2;

	imA = resizeX(image1, nx1, ny1, nz1, dx1, nx2, dx2);

	imB = resizeY(imA, nx2, ny1, nz1, dy1, ny2, dy2);
	free(imA);
   
	image2=resizeZ(imB, nx2, ny2, nz1, dz1, nz2, dz2);
	free(imB);

	return(image2);
}

short *resizeXYZ(short *image1,  DIM dim1, DIM dim2)
{
	float *imA, *imB;
	short *image2;

	imA = resizeX(image1, dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim2.nx, dim2.dx);

	imB = resizeY(imA, dim2.nx, dim1.ny, dim1.nz, dim1.dy, dim2.ny, dim2.dy);
	free(imA);
   
	image2=resizeZ(imB, dim2.nx, dim2.ny, dim1.nz, dim1.dz, dim2.nz, dim2.dz);
	free(imB);

	return(image2);
}

short *resizeZ(float *image1, int nx1, int ny1, int nz1, float dz1, int nz2, float dz2)
{
   int n;
   int np1;
   int k0;
   float d,dl;
   float sd=0.0;
   float *h;
   float *x;
   float zc1,zc2;
   short *image2;

   np1=nx1*ny1;

   zc1 = dz1*(nz1-1.0)/2.0;
   zc2 = dz2*(nz2-1.0)/2.0;

   image2=(short *)calloc(np1*nz2,sizeof(short));
   if(image2==NULL) return(NULL);

   if(dz1 < dz2)
   {
      sd=(float)sqrt( (0.5/log(2.0))*( dz2*dz2 - dz1*dz1 ) );
      sd /= dz1;
   }
   else sd = 0.0;

   h = gaussian_kernel(sd,&n);

   x=(float *)calloc(nz1,sizeof(float));

   for(int j=0;j<ny1;j++)
   for(int i=0;i<nx1;i++)
   {
      for(int l=0;l<nz1;l++)
         x[l]=image1[np1*l + nx1*j +i];

      for(int k=0;k<nz2;k++)
      {
         d = (k*dz2 - zc2 + zc1)/dz1;
         k0=(int)d;
         dl=d - k0;
         
         if(dl==0)
            image2[np1*k+nx1*j+i]=(short)(conv_pnt_sk(x,nz1,h,n,k0)+0.5);
         else
            image2[np1*k+nx1*j+i]=(short) ( (1.0-dl)*conv_pnt_sk(x,nz1,h,n,k0) +
            dl*conv_pnt_sk(x,nz1,h,n,k0+1) + 0.5 );
      }
   }

   free(h);
   free(x);

   return(image2);
}

float *resizeXYZ(float *image1, 
	int nx1, int ny1, int nz1, 
	float dx1, float dy1, float dz1,
	int nx2, int ny2, int nz2, 
	float dx2, float dy2, float dz2)
{
	float *imA, *imB;
	float *image2;

	imA = resizeX(image1, nx1, ny1, nz1, dx1, nx2, dx2);

	imB = resizeY(imA, nx2, ny1, nz1, dy1, ny2, dy2);
	free(imA);
   
	image2=resizeZ(imB, nx2, ny2, nz1, dz1, &nz2, &dz2);
	free(imB);

	return(image2);
}

float *resizeY(float *image1, int nx1, int ny1, int nz1, float dy1, int ny2, float dy2)
{
   int n;
   int np1, np2;
   int j0;
   float d,dl;
   float sd=0.0;
   float *h;
   float *image2;
   float *x;
   float yc1,yc2;

   np1=nx1*ny1;
   np2=nx1*ny2;

   yc1 = dy1*(ny1-1.0)/2.0;
   yc2 = dy2*(ny2-1.0)/2.0;

   image2=(float *)calloc(np2*nz1,sizeof(float));
   if(image2==NULL) return(NULL);

   if(dy1 < dy2)
   {
      sd=(float)sqrt( (0.5/log(2.0))*( dy2*dy2 - dy1*dy1 ) );
      sd /= dy1;
   }
   else sd=0.0;

   h = gaussian_kernel(sd,&n);

   x=(float *)calloc(ny1,sizeof(float));

   for(int k=0;k<nz1;k++)
   for(int i=0;i<nx1;i++)
   {
      for(int l=0;l<ny1;l++)
         x[l]=image1[np1*k + nx1*l +i];

      for(int j=0;j<ny2;j++)
      {
         d = (j*dy2 - yc2 + yc1)/dy1;
         j0=(int)d;
         dl=d - j0;

         if(dl==0.0)
            image2[np2*k+nx1*j+i]=conv_pnt_sk(x,ny1,h,n,j0);
         else
            image2[np2*k+nx1*j+i]=(1.0-dl)*conv_pnt_sk(x,ny1,h,n,j0) + dl*conv_pnt_sk(x,ny1,h,n,j0+1);
      }
   }

   free(h);
   free(x);

   return(image2);
}

float *resizeZ(float *image1, int nx1, int ny1, int nz1, float dz1, int *nz2, float *dz2)
{
   int n;
   int np1;
   int k0;
   float d,dl;
   float sd=0.0;
   float *h;
   float *x;
   float zc1,zc2;
   float *image2;

   np1=nx1*ny1;

   zc1 = dz1*(nz1-1.0)/2.0;
   zc2 = (*dz2)*( (*nz2) -1.0)/2.0;

   image2=(float *)calloc(np1*(*nz2) ,sizeof(float));
   if(image2==NULL) return(NULL);

   if( dz1 < (*dz2) )
   {
      sd=(float)sqrt( (0.5/log(2.0))*( (*dz2)*(*dz2) - dz1*dz1 ) );
      sd /= dz1;
   } else sd=0.0;

   h = gaussian_kernel(sd,&n);

   x=(float *)calloc(nz1,sizeof(float));

   for(int j=0;j<ny1;j++)
   for(int i=0;i<nx1;i++)
   {
      for(int l=0;l<nz1;l++)
         x[l]=image1[np1*l + nx1*j +i];

      for(int k=0;k<(*nz2);k++)
      {
         d = (k*(*dz2) - zc2 + zc1)/dz1;
         k0=(int)d;
         dl=d - k0;

         if(dl==0)
            image2[np1*k+nx1*j+i]=conv_pnt_sk(x,nz1,h,n,k0);
         else
            image2[np1*k+nx1*j+i]= (1.0-dl)*conv_pnt_sk(x,nz1,h,n,k0) + dl*conv_pnt_sk(x,nz1,h,n,k0+1);
      }
   }

   free(h);
   free(x);

   return(image2);
}
