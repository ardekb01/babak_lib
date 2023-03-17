#define _utils

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <spm_analyze.h>
#include <niftiimage.h>
#include <babak_lib.h>
#include <smooth.h>

void art_to_fsl(float *Mart, float *Mfsl, DIM sub_dim, DIM trg_dim);
void fsl_to_art(float *Mfsl, float *Mart, DIM sub_dim, DIM trg_dim);
int ccsize(short *im, int nv);
void checkDimension(int N, char **imagefile, int nx, int ny, int nz);
void affineLSE(char *msk, int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
void affineLSE(short *msk, int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
float *affineLSE(char *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp);
void affineLSE(char *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp, float *T);
void extractArray(short *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
void extractArray(float *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
void extractArray(short *im, int nx, int ny, int nx0, int ny0, int Lx, int Ly, float *array);
void extractArray(short *im, int nx, int i, int j, int L, float *array);
void getfilename(char *filename, const char *path);
void centerOfMass(short *im, int nx, int ny, int nz, float dx, float dy, float dz, float *CM);
void printMatrix(float *mat, int n, int p, const char *s, FILE *fp);
void printMatrix(int *mat, int n, int p, const char *s, FILE *fp);
void printMatrix(double *mat, int n, int p, const char *s, FILE *fp);
void get_temp_filename(char *filename);
void mask_and_save(const char *inputfile, const char *outputfile, short *mask, short *masked_image, int nbv, float FWHM);
void mask_and_save_nii(const char *inputfile, const char *outputfile, short *mask, short *masked_image, int nbv, float FWHM);
void read_transpose_save(char *inputfile, char *outputfile, int nr, int v);
short *readMask(const char *filename, int *nx, int *ny, int *nz);
float *readDataMatrix(char **imageList, int n, int p, short *mask);
void sobel_edge(short *in, float *out, int nx, int ny);

//////////////////////////////////////////////////////////////////////////////////

void art_to_fsl(float *Mart, float *Mfsl, DIM sub_dim, DIM trg_dim)
{
   float Tsub[16], Ttrg[16];
   float *inv_Tsub;

   Tsub[0]=1.0;  Tsub[1]=0.0;  Tsub[2]=0.0;  Tsub[3]=(sub_dim.nx-1.0)*sub_dim.dx/2.0;
   Tsub[4]=0.0;  Tsub[5]=1.0;  Tsub[6]=0.0;  Tsub[7]=(sub_dim.ny-1.0)*sub_dim.dy/2.0;
   Tsub[8]=0.0;  Tsub[9]=0.0;  Tsub[10]=1.0; Tsub[11]=(sub_dim.nz-1.0)*sub_dim.dz/2.0;
   Tsub[12]=0.0; Tsub[13]=0.0; Tsub[14]=0.0; Tsub[15]=1.0;

   Ttrg[0]=1.0;  Ttrg[1]=0.0;  Ttrg[2]=0.0;  Ttrg[3]=(trg_dim.nx-1.0)*trg_dim.dx/2.0;
   Ttrg[4]=0.0;  Ttrg[5]=1.0;  Ttrg[6]=0.0;  Ttrg[7]=(trg_dim.ny-1.0)*trg_dim.dy/2.0;
   Ttrg[8]=0.0;  Ttrg[9]=0.0;  Ttrg[10]=1.0; Ttrg[11]=(trg_dim.nz-1.0)*trg_dim.dz/2.0;
   Ttrg[12]=0.0; Ttrg[13]=0.0; Ttrg[14]=0.0; Ttrg[15]=1.0;

   inv_Tsub = inv4(Tsub);

   multi(Ttrg,4,4,Mart,4,4,Mfsl);
   multi(Mfsl,4,4,inv_Tsub,4,4, Mfsl);

   free(inv_Tsub);
}

void fsl_to_art(float *Mfsl, float *Mart, DIM sub_dim, DIM trg_dim)
{
   float Tsub[16], Ttrg[16];
   float *inv_Ttrg;

   Tsub[0]=1.0;  Tsub[1]=0.0;  Tsub[2]=0.0;  Tsub[3]=(sub_dim.nx-1.0)*sub_dim.dx/2.0;
   Tsub[4]=0.0;  Tsub[5]=1.0;  Tsub[6]=0.0;  Tsub[7]=(sub_dim.ny-1.0)*sub_dim.dy/2.0;
   Tsub[8]=0.0;  Tsub[9]=0.0;  Tsub[10]=1.0; Tsub[11]=(sub_dim.nz-1.0)*sub_dim.dz/2.0;
   Tsub[12]=0.0; Tsub[13]=0.0; Tsub[14]=0.0; Tsub[15]=1.0;

   Ttrg[0]=1.0;  Ttrg[1]=0.0;  Ttrg[2]=0.0;  Ttrg[3]=(trg_dim.nx-1.0)*trg_dim.dx/2.0;
   Ttrg[4]=0.0;  Ttrg[5]=1.0;  Ttrg[6]=0.0;  Ttrg[7]=(trg_dim.ny-1.0)*trg_dim.dy/2.0;
   Ttrg[8]=0.0;  Ttrg[9]=0.0;  Ttrg[10]=1.0; Ttrg[11]=(trg_dim.nz-1.0)*trg_dim.dz/2.0;
   Ttrg[12]=0.0; Ttrg[13]=0.0; Ttrg[14]=0.0; Ttrg[15]=1.0;

   inv_Ttrg = inv4(Ttrg);

   multi(inv_Ttrg,4,4,Mfsl,4,4,Mart);
   multi(Mart,4,4,Tsub,4,4, Mart);

   free(inv_Ttrg);
}

//////////////////////////////////////////////////////////////////////////////////

// Set the nxn matrix A equal to the identity matrix
void set_to_I( float *A, int n)
{
   for(int i=0; i<n*n; i++) 
   {
      A[i]=0.0;
   }

   for(int i=0; i<n; i++)
   {
      A[n*i + i] = 1.0;
   }
}

//////////////////////////////////////////////////////////////////////////////////

int ccsize(short *im, int nv)
{
   int count=0;

   if(im==NULL) return(count);

   for(int v=0; v<nv; v++)
   {
      if(im[v]>0) 
      {
         count++;
      }
   }

   return(count);
}

//////////////////////////////////////////////////////////////////////////////////
void copyarray(float *source, float *destination, int size)
{
   if(source == NULL || destination==NULL || size<=0) 
   {
      return;
   }

   for(int i=0; i<size; i++)
   {
      destination[i] = source[i];
   }

   return;
}

void copyarray(short *source, short *destination, int size)
{
   if(source == NULL || destination==NULL || size<=0) 
   {
      return;
   }

   for(int i=0; i<size; i++)
   {
      destination[i] = source[i];
   }

   return;
}

void copyarray(short *source, char *destination, int size)
{
   if(source == NULL || destination==NULL || size<=0) 
   {
      return;
   }

   for(int i=0; i<size; i++)
   {
      destination[i] = (char)(source[i]);
   }

   return;
}

void copyarray(char *source, short *destination, int size)
{
   if(source == NULL || destination==NULL || size<=0) 
   {
      return;
   }

   for(int i=0; i<size; i++)
   {
      destination[i] = (short)(source[i]);
   }

   return;
}

void zeroarray(float *y, int size)
{
   if(y == NULL || size<=0 ) 
   {
      return;
   }

   for(int i=0; i<size; i++)
   {
      y[i] = 0.0;
   }

   return;
}
//////////////////////////////////////////////////////////////////////////////////

float diceindex(short *setA, short *setB, int n)
{
   int sizeA=0, sizeB=0, sizeAandB=0;

   if(n<=0) 
   {
      printf("\n\nWarning: A non-positive array dimension passed to the diceindex() function.\n\n");
      return(0.0);
   }

   if( setA==NULL || setB==NULL) 
   {
      printf("\n\nWarning: A NULL array passed to the diceindex() function.\n\n");
      return(0.0);
   }

   for(int i=0; i<n; i++)
   {
      if( setA[i] > 0 )
         sizeA++;

      if( setB[i] > 0 )
         sizeB++;

      if( setA[i] > 0 && setB[i] > 0)
         sizeAandB++;
   }

   if( (sizeA + sizeB) == 0 )
      return(0.0);

   return( 2.0*sizeAandB/(sizeA+sizeB) );
}
//////////////////////////////////////////////////////////////////////////////////

void saveMatrix(float *A, int n, int m, char *filename)
{
   FILE *fp;

   fp = fopen(filename, "w");
   if(fp==NULL) file_open_error(filename);

   fwrite(&n, sizeof(int), 1, fp);
   fwrite(&m, sizeof(int), 1, fp);
   fwrite(A, sizeof(float), n*m, fp);

   fclose(fp);
}

void saveMatrix(short *A, int n, int m, char *filename)
{
   FILE *fp;

   fp = fopen(filename, "w");
   if(fp==NULL) file_open_error(filename);

   fwrite(&n, sizeof(int), 1, fp);
   fwrite(&m, sizeof(int), 1, fp);
   fwrite(A, sizeof(short), n*m, fp);

   fclose(fp);
}

//////////////////////////////////////////////////////////////////////////////////
float *readMatrix(int *n, int *m, char *filename)
{
   FILE *fp;
   float *A;

   fp = fopen(filename, "r");
   if(fp==NULL) file_open_error(filename);

   if( fread(n, sizeof(int), 1, fp) != 1)
   {
      printf("Failed to read variable \"n\".\n");
   }

   if( fread(m, sizeof(int), 1, fp) != 1)
   {
      printf("Failed to read variable \"m\".\n");
   }

   A = (float *)calloc( (*n)*(*m), sizeof(float));

   if( fread(A, sizeof(float), (*n)*(*m), fp) != (size_t)((*n)*(*m)) ) 
   {
      printf("Failed to read variable \"A\".\n");
   }

   fclose(fp);

   return(A);
}
//////////////////////////////////////////////////////////////////////////////////

// returns 1 if all images have the same dimensions nx, ny, and nz, 0 otherwise
void checkDimension(int N, char **imagefile, int nx, int ny, int nz)
{
	struct dsr analyzehdr;
	char hdrfile[1024];
	char imgfile[1024];
	int nx0, ny0, nz0;
	float dx,dy,dz;
	short dataType;

	if(N==0) return;

	for(int i=0; i<N; i++)
	{
		nx0=ny0=nz0=0;
		get_analyze_file_names(imagefile[i], hdrfile, imgfile);
		read_analyze_hdr(&analyzehdr, hdrfile);
		setDimensions(analyzehdr, &nx0, &ny0, &nz0, &dx, &dy, &dz, &dataType,0);

		if(nx != nx0 || ny != ny0 || nz != nz0) 
		{
			printf("\n\nImage %d: %s",i+1,imagefile[i]);
			printf("\n\tMatrix size = %d x %d x %d", nx0,ny0,nz0);
			printf("\n\nAll input images must be of size: %d x %d x %d, aborting ...\n\n",nx,ny,nz);
			exit(0);
		}
	}
}

void checkDimension_nifti(int N, char **imagefile, int nx, int ny, int nz)
{
   nifti_1_header hdr;

	if(N==0) return;

	for(int i=0; i<N; i++)
	{
		hdr = read_NIFTI_hdr(imagefile[i]);

		if(nx != hdr.dim[1] || ny != hdr.dim[2] || nz != hdr.dim[3]) 
		{
			printf("\n\nImage %d: %s",i+1,imagefile[i]);
			printf("\n\tMatrix size = %d x %d x %d", hdr.dim[1],hdr.dim[2],hdr.dim[3]);
			printf("\n\nAll input images must be of size: %d x %d x %d, aborting ...\n\n",nx,ny,nz);
			exit(0);
		}
	}
}

void checkDimension(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
	struct dsr analyzehdr;
	char hdrfile[1024];
	char imgfile[1024];
	int nx0, ny0, nz0;
	float dx0,dy0,dz0;
	short dataType;

	if(N==0) return;

	get_analyze_file_names(imagefile[0], hdrfile, imgfile);
	read_analyze_hdr(&analyzehdr, hdrfile);
	setDimensions(analyzehdr, nx, ny, nz, dx, dy, dz, &dataType,0);

	for(int i=0; i<N; i++)
	{
		nx0=ny0=nz0=0;
		get_analyze_file_names(imagefile[i], hdrfile, imgfile);
		read_analyze_hdr(&analyzehdr, hdrfile);
		setDimensions(analyzehdr, &nx0, &ny0, &nz0, &dx0, &dy0, &dz0, &dataType,0);

		if( *nx != nx0 ||  *ny != ny0 ||  *nz != nz0) 
		{
			printf("\n\nImage %d: %s",i+1,imagefile[i]);
			printf("\n\tMatrix size = %d x %d x %d", nx0,ny0,nz0);
			printf("\n\nAborting ...\n\n");
			exit(0);
		}
	}
}

void affineLSE(char *msk, int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
	int np;
	int q,N;
	float xc, yc, zc;
	float rx,ry,rz;
	float sx,sy,sz;
  	float  x,y,z;   
	float *AFF, *invAFF;	// affine transform
	float *invT;

	double Mrx, Mry, Mrz;
	double Msx, Msy, Msz;
	double SR[9], RR[9];
	double *invRR;
	double A[9],B[3];

	AFF = (float *)calloc(16,sizeof(float));

	invT = inv4(T);

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

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = invT[0]*x + invT[1]*y + invT[2]*z  + invT[3];
			sy = invT[4]*x + invT[5]*y + invT[6]*z  + invT[7];
			sz = invT[8]*x + invT[9]*y + invT[10]*z + invT[11];

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

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = invT[0]*x + invT[1]*y + invT[2]*z  + invT[3];
			sy = invT[4]*x + invT[5]*y + invT[6]*z  + invT[7];
			sz = invT[8]*x + invT[9]*y + invT[10]*z + invT[11];

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
	// This portion of code adjusts the warp field components according to Eq. (4)
	// in the technical notes.
	///////////////////////////////////////////////////////////////////////////////
	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		// find r at point (i,j,k)
		rx = i*dx - xc;
		ry = j*dy - yc;
		rz = k*dz - zc;

		x = rx + Xwarp[q];
		y = ry + Ywarp[q];
		z = rz + Zwarp[q];

		// find s the image of r in the object space
		sx = invT[0]*x + invT[1]*y + invT[2]*z  + invT[3];
		sy = invT[4]*x + invT[5]*y + invT[6]*z  + invT[7];
		sz = invT[8]*x + invT[9]*y + invT[10]*z + invT[11];

		// adjust warp field components
		Xwarp[q] = invAFF[0]*sx + invAFF[1]*sy + invAFF[2]*sz + invAFF[3] - rx;
		Ywarp[q] = invAFF[4]*sx + invAFF[5]*sy + invAFF[6]*sz + invAFF[7] - ry;
		Zwarp[q] = invAFF[8]*sx + invAFF[9]*sy + invAFF[10]*sz + invAFF[11] - rz;
	}
	///////////////////////////////////////////////////////////////////////////////

	// replace T with invAFF
	for(int i=0; i<16; i++) T[i] = invAFF[i];

	delete invAFF;
}

void affineLSE(short *msk, int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
	int np;
	int q,N;
	float xc, yc, zc;
	float rx,ry,rz;
	float sx,sy,sz;
  	float  x,y,z;   
	float *AFF, *invAFF;	// affine transform
	float *invT;

	double Mrx, Mry, Mrz;
	double Msx, Msy, Msz;
	double SR[9], RR[9];
	double *invRR;
	double A[9],B[3];

	AFF = (float *)calloc(16,sizeof(float));

	invT = inv4(T);

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

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = invT[0]*x + invT[1]*y + invT[2]*z  + invT[3];
			sy = invT[4]*x + invT[5]*y + invT[6]*z  + invT[7];
			sz = invT[8]*x + invT[9]*y + invT[10]*z + invT[11];

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

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = invT[0]*x + invT[1]*y + invT[2]*z  + invT[3];
			sy = invT[4]*x + invT[5]*y + invT[6]*z  + invT[7];
			sz = invT[8]*x + invT[9]*y + invT[10]*z + invT[11];

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
	// This portion of code adjusts the warp field components according to Eq. (4)
	// in the technical notes.
	///////////////////////////////////////////////////////////////////////////////
	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		// find r at point (i,j,k)
		rx = i*dx - xc;
		ry = j*dy - yc;
		rz = k*dz - zc;

		x = rx + Xwarp[q];
		y = ry + Ywarp[q];
		z = rz + Zwarp[q];

		// find s the image of r in the object space
		sx = invT[0]*x + invT[1]*y + invT[2]*z  + invT[3];
		sy = invT[4]*x + invT[5]*y + invT[6]*z  + invT[7];
		sz = invT[8]*x + invT[9]*y + invT[10]*z + invT[11];

		// adjust warp field components
		Xwarp[q] = invAFF[0]*sx + invAFF[1]*sy + invAFF[2]*sz + invAFF[3] - rx;
		Ywarp[q] = invAFF[4]*sx + invAFF[5]*sy + invAFF[6]*sz + invAFF[7] - ry;
		Zwarp[q] = invAFF[8]*sx + invAFF[9]*sy + invAFF[10]*sz + invAFF[11] - rz;
	}
	///////////////////////////////////////////////////////////////////////////////

   // replace T with invAFF
   for(int i=0; i<16; i++) T[i] = invAFF[i];

   free(invT);
   delete invAFF;
}

float *affineLSE(char *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp)
{
	int q,N;
	float xc, yc;
	float rx,ry;
	float sx,sy;
  	//float  x,y;   
	float *T;

	double Mrx, Mry;
	double Msx, Msy;
	double SR[4], RR[4];
	double *invRR;
	double A[4],B[2];

	T = (float *)calloc(9,sizeof(float));

   	xc=dx*(nx-1)/2.0;     /* +---+---+ */
	yc=dy*(ny-1)/2.0;

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes averages of the s and r vectors defined in
	// the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	Mrx=Mry=0.0;
	Msx=Msy=0.0;

	N = 0;
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = j*nx + i;

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;

			sx = rx + Xwarp[q];
			sy = ry + Ywarp[q];

			Mrx += rx; Mry += ry;
			Msx += sx; Msy += sy;

			N++;
		}
	}

	if(N!=0)
	{
		Mrx /= N; Mry /= N;
		Msx /= N; Msy /= N;
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes the two 2x2 matrix in Eq. (2) of the
	// tech. notes.
	/////////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<4; i++) SR[i]=RR[i]=0.0;

	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = j*nx + i;

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;

			sx = rx + Xwarp[q];
			sy = ry + Ywarp[q];

			rx -= Mrx; ry -= Mry;
			sx -= Msx; sy -= Msy;

			SR[0]+=sx*rx; SR[1]+=sx*ry; 
			SR[2]+=sy*rx; SR[3]+=sy*ry; 

			RR[0]+=rx*rx; RR[1]+=rx*ry; 
			RR[2]+=ry*rx; RR[3]+=ry*ry;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// estimate A according to Eq. (2) of the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	invRR = inv2(RR);
	multi(SR,2,2,invRR,2,2,A);
	free(invRR);
	/////////////////////////////////////////////////////////////////////////////////

	// estimate B according to Eq. (1) of the technical notes
	B[0] = Msx - A[0]*Mrx - A[1]*Mry;
	B[1] = Msy - A[2]*Mrx - A[3]*Mry;

	// Eq. (3) of tech. notes
	T[0]=(float)A[0]; T[1]=(float)A[1];  T[2]=(float)B[0];
	T[3]=(float)A[2]; T[4]=(float)A[3];  T[5]=(float)B[1]; 
	T[6]=0.0; T[7]=0.0; T[8]=1.0;

	return(T);
}
//////////////////////////////////////////////////////////////////////////////////////////////////

void affineLSE(char *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp, float *T)
{
	//int np;
	int q,N;
	float xc, yc;
	float rx,ry;
	float sx,sy;
  	//float  x,y;   

	double Mrx, Mry;
	double Msx, Msy;
	double SR[4], RR[4];
	double *invRR;
	double A[4],B[2];

	//np = nx*ny;

   	xc=dx*(nx-1)/2.0;     /* +---+---+ */
	yc=dy*(ny-1)/2.0;

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes averages of the s and r vectors defined in
	// the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	Mrx=Mry=0.0;
	Msx=Msy=0.0;

	N = 0;
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = j*nx + i;

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;

			sx = rx + Xwarp[q];
			sy = ry + Ywarp[q];

			Mrx += rx; Mry += ry;
			Msx += sx; Msy += sy;

			N++;
		}
	}

	if(N!=0)
	{
		Mrx /= N; Mry /= N;
		Msx /= N; Msy /= N;
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes the two 2x2 matrix in Eq. (2) of the
	// tech. notes.
	/////////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<4; i++) SR[i]=RR[i]=0.0;

	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = j*nx + i;

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;

			sx = rx + Xwarp[q];
			sy = ry + Ywarp[q];

			rx -= Mrx; ry -= Mry;
			sx -= Msx; sy -= Msy;

			SR[0]+=sx*rx; SR[1]+=sx*ry; 
			SR[2]+=sy*rx; SR[3]+=sy*ry; 

			RR[0]+=rx*rx; RR[1]+=rx*ry; 
			RR[2]+=ry*rx; RR[3]+=ry*ry;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// estimate A according to Eq. (2) of the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	invRR = inv2(RR);
	multi(SR,2,2,invRR,2,2,A);
	free(invRR);
	/////////////////////////////////////////////////////////////////////////////////

	// estimate B according to Eq. (1) of the technical notes
	B[0] = Msx - A[0]*Mrx - A[1]*Mry;
	B[1] = Msy - A[2]*Mrx - A[3]*Mry;

	// Eq. (3) of tech. notes
	T[0]=(float)A[0]; T[1]=(float)A[1];  T[2]=(float)B[0];
	T[3]=(float)A[2]; T[4]=(float)A[3];  T[5]=(float)B[1]; 
	T[6]=0.0; T[7]=0.0; T[8]=1.0;
}

void affineLSE(short *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp, float *T)
{
	//int np;
	int q,N;
	float xc, yc;
	float rx,ry;
	float sx,sy;
  	//float  x,y;   

	double Mrx, Mry;
	double Msx, Msy;
	double SR[4], RR[4];
	double *invRR;
	double A[4],B[2];

	//np = nx*ny;

   	xc=dx*(nx-1)/2.0;     /* +---+---+ */
	yc=dy*(ny-1)/2.0;

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes averages of the s and r vectors defined in
	// the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	Mrx=Mry=0.0;
	Msx=Msy=0.0;

	N = 0;
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = j*nx + i;

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;

			sx = rx + Xwarp[q];
			sy = ry + Ywarp[q];

			Mrx += rx; Mry += ry;
			Msx += sx; Msy += sy;

			N++;
		}
	}

	if(N!=0)
	{
		Mrx /= N; Mry /= N;
		Msx /= N; Msy /= N;
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes the two 2x2 matrix in Eq. (2) of the
	// tech. notes.
	/////////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<4; i++) SR[i]=RR[i]=0.0;

	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = j*nx + i;

		if(msk[q])
		{
			rx = i*dx - xc;
			ry = j*dy - yc;

			sx = rx + Xwarp[q];
			sy = ry + Ywarp[q];

			rx -= Mrx; ry -= Mry;
			sx -= Msx; sy -= Msy;

			SR[0]+=sx*rx; SR[1]+=sx*ry; 
			SR[2]+=sy*rx; SR[3]+=sy*ry; 

			RR[0]+=rx*rx; RR[1]+=rx*ry; 
			RR[2]+=ry*rx; RR[3]+=ry*ry;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// estimate A according to Eq. (2) of the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	invRR = inv2(RR);
	multi(SR,2,2,invRR,2,2,A);
	free(invRR);
	/////////////////////////////////////////////////////////////////////////////////

	// estimate B according to Eq. (1) of the technical notes
	B[0] = Msx - A[0]*Mrx - A[1]*Mry;
	B[1] = Msy - A[2]*Mrx - A[3]*Mry;

	// Eq. (3) of tech. notes
	T[0]=(float)A[0]; T[1]=(float)A[1];  T[2]=(float)B[0];
	T[3]=(float)A[2]; T[4]=(float)A[3];  T[5]=(float)B[1]; 
	T[6]=0.0; T[7]=0.0; T[8]=1.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void extractArray(short *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array)
{
   int q=0;
   int np;
   int znp, ynx;

   np=nx*ny;
   for(int z=nz0-Lz; z<=(nz0+Lz); z++)
   {
      znp = z*np;
      for(int y=ny0-Ly; y<=(ny0+Ly); y++)
      {
         ynx = y*nx;
         for(int x=nx0-Lx; x<=(nx0+Lx); x++)
         {
            if(y<0 || y>=ny || x<0 || x>=nx || z<0 || z>=nz) 
               array[q++]=0.0;
            else 
               array[q++]=im[znp + ynx + x]; 
         }
      }
   }
}

void extractArray(short *im, int nx, int ny, int nz, int np, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array)
{
   int q=0;
   int znp, ynx;

   for(int z=nz0-Lz; z<=(nz0+Lz); z++)
   {
      znp = z*np;
      for(int y=ny0-Ly; y<=(ny0+Ly); y++)
      {
         ynx = y*nx;
         for(int x=nx0-Lx; x<=(nx0+Lx); x++)
         {
            if(y<0 || y>=ny || x<0 || x>=nx || z<0 || z>=nz) 
               array[q++]=0.0;
            else 
               array[q++]=im[znp + ynx + x]; 
         }
      }
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void extractArray(float *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array)
{
   int q=0;
   int np;
   int znp, ynx;

   np=nx*ny;
   for(int z=nz0-Lz; z<=(nz0+Lz); z++)
   {
      znp = z*np;
      for(int y=ny0-Ly; y<=(ny0+Ly); y++)
      {
         ynx = y*nx;
         for(int x=nx0-Lx; x<=(nx0+Lx); x++)
         {
            if(y<0 || y>=ny || x<0 || x>=nx || z<0 || z>=nz) 
               array[q++]=0.0;
            else 
               array[q++]=im[znp + ynx + x]; 
         }
      }
   }
}

/*
void extractArray(float *im, int nx, int ny, int nz, int np, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array)
{
   int q=0;
   int znp, ynx;

   for(int z=nz0-Lz; z<=(nz0+Lz); z++)
   {
      znp = z*np;
      for(int y=ny0-Ly; y<=(ny0+Ly); y++)
      {
         ynx = y*nx;
         for(int x=nx0-Lx; x<=(nx0+Lx); x++)
         {
            if(y<0 || y>=ny || x<0 || x>=nx || z<0 || z>=nz) 
               array[q++]=0.0;
            else 
               array[q++]=im[znp + ynx + x]; 
         }
      }
   }
}
*/

// highly optimized version of the above
void extractArray(float *im, int nx, int ny, int nz, int np, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array)
{
   int q=0;
   int znp, ynx;
   int zi, zf;
   int yi, yf;
   int xi, xf;

   zi = nz0-Lz;
   zf = nz0+Lz;
   yi = ny0-Ly;
   yf = ny0+Ly;
   xi = nx0-Lx;
   xf = nx0+Lx;

   np=nx*ny;
   for(int z=zi; z<=zf; z++)
   {
      if(z<0 || z>=nz)
      {
         for(int y=yi; y<=yf; y++)
         for(int x=xi; x<=xf; x++)
               array[q++]=0.0;
         continue;
      }

      znp = z*np;
      for(int y=yi; y<=yf; y++)
      {
         if(y<0 || y>=ny) 
         {
            for(int x=xi; x<=xf; x++)
               array[q++]=0.0;

            continue;
         }

         ynx = y*nx;
         for(int x=xi; x<=xf; x++)
         {
            if(x<0 || x>=nx) 
               array[q++]=0.0;
            else 
               array[q++]=im[znp + ynx + x]; 
         }
      }
   }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void extractArray(short *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, short *array)
{
   int q=0;
   int np;
   int znp, ynx;

   np=nx*ny;
   for(int z=nz0-Lz; z<=(nz0+Lz); z++)
   {
      znp = z*np;
      for(int y=ny0-Ly; y<=(ny0+Ly); y++)
      {
         ynx = y*nx;
         for(int x=nx0-Lx; x<=(nx0+Lx); x++)
         {
            if(z<0 || z>=nz || y<0 || y>=ny || x<0 || x>=nx) 
               array[q++]=0;
            else 
               array[q++]=im[znp + ynx + x]; 
         }
      }
   }
}

// rv will be one if a zero value is encountered
// highly optimized version of the above
int extractArray(short *im, int nx, int ny, int nz, int np, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, short *array)
{
   int rv=0; // return value
   short voxel_val;

   int q=0;
   int znp, ynx;
   int zi, zf;
   int yi, yf;
   int xi, xf;

   zi = nz0-Lz;
   zf = nz0+Lz;
   yi = ny0-Ly;
   yf = ny0+Ly;
   xi = nx0-Lx;
   xf = nx0+Lx;

   np=nx*ny;
   for(int z=zi; z<=zf; z++)
   {
      if(z<0 || z>=nz)
      {
         for(int y=yi; y<=yf; y++)
         for(int x=xi; x<=xf; x++)
               array[q++]=0;

         rv=1;
         continue;
      }

      znp = z*np;
      for(int y=yi; y<=yf; y++)
      {
         if(y<0 || y>=ny) 
         {
            for(int x=xi; x<=xf; x++)
               array[q++]=0;

            rv=1;
            continue;
         }

         ynx = y*nx;
         for(int x=xi; x<=xf; x++)
         {
            if(x<0 || x>=nx) 
            {
               array[q++]=0;
               rv=1;
            }
            else 
            {
               voxel_val = im[znp + ynx + x]; 
               array[q++]=voxel_val; 
               if(voxel_val==0) rv=1;
            }
         }
      }
   }

   return(rv);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// nx, ny, nz: im matrix dimensions
// x0, y0, z0: template center coordinates
// n: number of voxels in the template
// x, y, z: template voxel coordinates wrt (x0,y0,z0)
void extractArray(short *im,int nx,int ny,int nz,int x0,int y0,int z0,short *x,short *y,short *z,int n,float *array)
{
	int i,j,k;
	int np;

	np = nx*ny;

	for(int v=0; v<n; v++)
	{
		i = x0 + x[v];
		j = y0 + y[v];
		k = z0 + z[v];

		if(i<0 || i>=nx || j<0 || j>=ny || k<0 || k>=nz) 
			array[v]=0.0;
		else 
			array[v]=im[k*np + j*nx + i]; 
	}
}

void extractArray(unsigned char *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array)
{
	int q=0;
	int np;

	np=nx*ny;

	for(int z=nz0-Lz; z<=(nz0+Lz); z++)
	for(int y=ny0-Ly; y<=(ny0+Ly); y++)
	for(int x=nx0-Lx; x<=(nx0+Lx); x++)
	{
		if(y<0 || y>=ny || x<0 || x>=nx || z<0 || z>=nz) 
			array[q++]=0.0;
		else 
			array[q++]=im[z*np + y*nx + x]; 
	}
}

void extractArray(short *im, int nx, int ny, int nx0, int ny0, int Lx, int Ly, float *array)
{
	int q=0;

	for(int y=ny0-Ly; y<=(ny0+Ly); y++)
	for(int x=nx0-Lx; x<=(nx0+Lx); x++)
	{
		if(y<0 || y>=ny || x<0 || x>=nx ) 
			array[q++]=0.0;
		else 
			array[q++]=im[y*nx + x]; 
	}
}

void extractArray(short *im, int nx, int i, int j, int L, float *array)
{
	int n;

	for(int k=-L; k<=L; k++)
	{
		n = i+k;
		if(n<0 || n>=nx ) 
			array[k+L]=0.0;
		else 
			array[k+L]=im[nx*j+n]; 
	}
}

#if 0
// This function extracts the "filename" from the full "path" string.
// Example: If path="/home/babak/testimages/test.img", then filname="test.img".
void getfilename(char *filename, const char *path)
{
   int i;
   int len;	// length of the path string
   int fnlen;	// length of the filename string
   int pos;	// position of the filename
   char dum[512];
	
   len=(int)strlen(path);

   i=len-1;
   while( i>=0 && path[i] != '/' )
   {
      i--;
   }
   pos=i+1;

   strcpy(dum,path+pos);

   fnlen=(int)strlen(dum);

   strncpy(filename,path+pos,fnlen);
   filename[fnlen]='\0';

   return;
}
#endif

// This function extracts the "filename" from the full "path" string.
// Example: If path="/home/babak/testimages/test.img", then filname="test.img".
void getfilename(char *filename, const char *path)
{
   int i;
   int len; // length of the path string
   int pos; // position of the filename
	
   len=(int)strlen(path);

   i=len-1;
   while( i>=0 && path[i] != '/' )
   {
      i--;
   }
   pos=i+1;

   strcpy(filename, path+pos);

   return;
}

void centerOfMass(short *im, int nx, int ny, int nz, float dx, float dy, float dz, float *CM)
{
	//double trace;
	double sumx, sumy, sumz;
	double M;		// total mass
	//float X,Y,Z;
	float xc,yc,zc;
	int np;
	short v;
	float dzk, dyj;

	xc = dx*(nx-1.0)/2.0;
	yc = dy*(ny-1.0)/2.0;
	zc = dz*(nz-1.0)/2.0;
 
	sumx = sumy = sumz = M = 0.0;
   
	np = nx*ny;

	for(int k=0; k<nz; k++)
    {
		dzk = dz*k;
		for(int j=0; j<ny; j++)
		{
			dyj = dy*j;
			for(int i=0; i<nx; i++)
			{
				v = im[i + j*nx + k*np];

				sumx += dx*i*v;
				sumy += dyj*v;
				sumz += dzk*v;
				M += v;
			}
		}
	}

	CM[0] = (float)(sumx/M);
	CM[1] = (float)(sumy/M);
	CM[2] = (float)(sumz/M);

	CM[0] -= xc;
	CM[1] -= yc;
	CM[2] -= zc;
}

void printMatrix(float *mat, int n, int p, const char *s, FILE *fp)
{
	if(fp!=NULL)
	{
		fprintf(fp,"\n# %s\n",s);
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<p; j++)
				fprintf(fp,"%0.6f\t",mat[p*i + j]);
			fprintf(fp,"\n");
		}
	}
	else
	{
		printf("%s\n",s);
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<p; j++)
				printf("%0.6f\t",mat[p*i + j]);
			printf("\n");
		}
	}
}

void printMatrix(int *mat, int n, int p, const char *s, FILE *fp)
{
	if(fp!=NULL)
	{
		fprintf(fp,"\n# %s\n",s);
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<p; j++)
				fprintf(fp,"%d\t",mat[p*i + j]);
			fprintf(fp,"\n");
		}
	}
	else
	{
		printf("%s\n",s);
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<p; j++)
				printf("%d\t",mat[p*i + j]);
			printf("\n");
		}
	}
}

void printMatrix(double *mat, int n, int p, const char *s, FILE *fp)
{
	if(fp!=NULL)
	{
		fprintf(fp,"\n# %s\n",s);
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<p; j++)
				fprintf(fp,"%lf\t",mat[p*i + j]);
			fprintf(fp,"\n");
		}
	}
	else
	{
		printf("%s\n",s);
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<p; j++)
				printf("%lf\t",mat[p*i + j]);
			printf("\n");
		}
	}
}

void get_temp_filename(char *filename)
{
	FILE *fp;
	char *homedir;
	static int n=0;

	n++;

	sprintf(filename,"DELETEME%d%d",getpid(),n);
	fp=fopen(filename,"w");
	if(fp != NULL)
	{
		fclose(fp);
		remove(filename);
		return;
	}

	sprintf(filename,"/usr/tmp/DELETEME%d%d",getpid(),n);
	fp=fopen(filename,"w");
	if(fp != NULL)
	{
		fclose(fp);
		remove(filename);
		return;
	}

	sprintf(filename,"/tmp/DELETEME%d%d",getpid(),n);
	fp=fopen(filename,"w");
	if(fp != NULL)
	{
		fclose(fp);
		remove(filename);
		return;
	}

	homedir=getenv("HOME");
	if(homedir != NULL)
	{
		sprintf(filename,"%s/DELETEME%d%d",homedir,getpid(),n);
	}
	else
	{
		sprintf(filename,"DELETEME%d%d",getpid(),n);
	}

	fp=fopen(filename,"w");
	if(fp != NULL)
	{
		fclose(fp);
		remove(filename);
		return;
	}
}

void mask_and_save(const char *inputfile, const char *outputfile, short *mask, short *masked_image, int nbv, float FWHM)
{
	FILE *fp;
	short *image;
	int nx,ny,nz,nv;
	float dx,dy,dz;
	int type;
	int index;
	float sdx,sdy,sdz;
	float SD;

	// Open the file in which the masked images will be saved in append mode.
	fp = fopen(outputfile,"a");
	if( fp == NULL) file_open_error(outputfile);

	// Read a new image.
	image=(short *)read_analyze_image(inputfile,&nx, &ny, &nz, &dx, &dy, &dz, &type,0);

	// If for some reason cannot read this image, just go to the next one.
	if(image==NULL) 
	{
		fclose(fp);
		return;
	}

	nv=nx*ny*nz;

	SD = (float)( 0.5*FWHM / sqrt(2.0*log(2.0) ));

	if(SD>0.0 && dx!=0.0 && dy!=0.0 && dz!=0.0)
	{
		sdx = SD/dx; sdy = SD/dy; sdz = SD/dz;
		float *temp;
		temp = smoothXYZ(image, nx, ny, nz, sdx, sdy, sdz);
		for(int i=0; i<nv; i++) image[i]=(short)(temp[i]+0.5);
		delete temp;
	}

	// The actual masking is done here.
	index=0;
	for(int v=0; v<nv; v++) if(mask[v]!=0) masked_image[index++] = image[v];

	// save the masked image.
	fwrite(masked_image,sizeof(short),nbv,fp); 

	delete image;

	fclose(fp);
}

void mask_and_save_nii(const char *inputfile, const char *outputfile, short *mask, short *masked_image, int nbv, float FWHM)
{
   NIFTIIMAGE im;
   short *image=NULL;

   FILE *fp;
   int nv;
   //int type;
   int index;
   float sdx,sdy,sdz;
   float dx,dy,dz;
   float SD;

   // Open the file in which the masked images will be saved in append mode.
   fp = fopen(outputfile,"a");
   if( fp == NULL) file_open_error(outputfile);

   // Read a new image.
   im.read(inputfile);
   image = (short *)im.getdata();

   if(image==NULL) 
   {
      fclose(fp);
      return;
   }

   nv=im.nv();
   dx = im.dx();
   dy = im.dy();
   dz = im.dz();

   SD = (float)( 0.5*FWHM / sqrt(2.0*log(2.0) ));

   if(SD>0.0 && dx!=0.0 && dy!=0.0 && dz!=0.0)
   {
      float *temp;

      sdx = SD/dx; sdy = SD/dy; sdz = SD/dz;
      temp = smoothXYZ(image, im.nx(), im.ny(), im.nz(), sdx, sdy, sdz);
      for(int i=0; i<nv; i++) image[i]=(short)(temp[i]+0.5);

      delete temp;
   }

   // The actual masking is done here.
   index=0;
   for(int v=0; v<nv; v++) if(mask[v]!=0) masked_image[index++] = image[v];

   // save the masked image.
   fwrite(masked_image,sizeof(short),nbv,fp); 

   fclose(fp);
}

void read_transpose_save(char *inputfile, char *outputfile, int nr, int v)
{
	short *X=NULL, *x=NULL; 
	int nc;	// number of columns of the input matrix
	struct stat fileinfo;   // structure into which information is placed about the file
	FILE *fpi,*fpo;
	int	index;

	// obtain information about the masked image series file.
	if( stat(inputfile, &fileinfo) == -1 )
	{
		printf("\n\nstat() failure. Cannot read %s.\n\n",inputfile);
		exit(0);
	}

	nc = fileinfo.st_size / (nr*sizeof(short));

	if(v) printf("Input matrix size = %d x %d\n",nr,nc);

	if( nc<=0 || nr<=0) 
	{
		printf("\n\nread_transpose_save(): Input matrix size is %d x %d, aborting ...\n\n",nr,nc);
		exit(0);
	}

	// allocated memory for masked images.
	X = (short *)calloc(nr*nc,sizeof(short));
	if(X==NULL) memory_allocation_error("read_transpose_save(): X");

	x = (short *)calloc(nr,sizeof(short));
	if(x==NULL) memory_allocation_error("read_transpose_save(): x");

	// Open the file from which the input matrix will be read.
	fpi = fopen(inputfile,"r");
	if( fpi == NULL) file_open_error(inputfile);

	index=fread(X,sizeof(short),nr*nc,fpi); 

	if(index!=(nr*nc))
	{
		printf("\n\nfread() failure. Cannot read %s, aborting ...\n",inputfile);
		exit(0);
	}

	fclose(fpi);

	// Open the file to which the output matrix will be written.
	fpo = fopen(outputfile,"w");
	if( fpo == NULL) file_open_error(outputfile);

	for(int j=0; j<nc; j++)
	{
		for(int i=0; i<nr; i++)
			x[i] = X[i*nc+j];

		fwrite(x,sizeof(short),nr,fpo); 
	}

	fclose(fpo);

	// free memory
	if(x!=NULL) delete x;	
	if(X!=NULL) delete X;	
}

short *readMask(const char *filename, int *nx, int *ny, int *nz)
{
	short *mask;
	float dx, dy, dz;
	int type;

	mask =(short *)read_analyze_image(filename, nx, ny, nz, &dx, &dy, &dz, &type, 0);

	if(mask==NULL)
	{
		printf("\nError: Reading image %s failed, aborting ...\n\n",filename);
		exit(0);
	}

	if(type != 4)
	{
		printf("\nError: Mask image is not of type \"short\", aborting ...\n\n");
		exit(0);
	}

	return(mask);
}

short *readMask_nifti(const char *filename, int *nx, int *ny, int *nz)
{
   nifti_1_header hdr;
   short *mask;

   mask = (short *)read_nifti_image(filename, &hdr);
  
   *nx = hdr.dim[1];
   *ny = hdr.dim[2];
   *nz = hdr.dim[3];

   if(mask==NULL)
   {
      printf("\nError: Reading image %s failed, aborting ...\n\n",filename);
      exit(0);
   }

   if(hdr.datatype != 4)
   {
      printf("\nError: Mask image is not of type \"short\", aborting ...\n\n");
      exit(0);
   }

   return(mask);
}

float *readDataMatrix(char **imageList, int n, int p, short *mask)
{
	short *im;
	int nx, ny, nz, nv, type;
	float dx,dy,dz;
	float *X;
	int j;

	X = (float *)calloc(n*p, sizeof(float));

	for(int i=0; i<n; i++)
	{
		j=0;

		im =(short *)read_analyze_image(imageList[i], &nx, &ny, &nz, &dx, &dy, &dz, &type, 0);
		nv = nx*ny*nz;

		if(im!=NULL)
		{
			for(int v=0; v<nv; v++)
			{
				if(mask[v]>0) 
					X[i*p + j++] = (float)im[v];
			}

			free(im);
		}
	}

	return(X);
}

float *readDataMatrix_nifti(char **imageList, int n, int p, short *mask)
{
	short *im;
	//int nx, ny, nz, nv, type;
	int nx, ny, nz, nv;
	//float dx,dy,dz;
	float *X;
	int j;
   nifti_1_header hdr;

	X = (float *)calloc(n*p, sizeof(float));

	for(int i=0; i<n; i++)
	{
		j=0;

		//im =(short *)read_analyze_image(imageList[i], &nx, &ny, &nz, &dx, &dy, &dz, &type, 0);
   		im = (short *)read_nifti_image(imageList[i], &hdr);
   		nx = hdr.dim[1];
		ny = hdr.dim[2];
  		nz = hdr.dim[3];
		nv = nx*ny*nz;

		if(im!=NULL)
		{
			for(int v=0; v<nv; v++)
			{
				if(mask[v]>0) 
					X[i*p + j++] = (float)im[v];
			}

			free(im);
		}
	}

	return(X);
}

short *readDataMatrixShort(char **imageList, int n, int p, short *mask)
{
	short *im;
	int nx, ny, nz, nv, type;
	float dx,dy,dz;
	short *X;
	int j;

	X = (short *)calloc(n*p, sizeof(short));

	for(int i=0; i<n; i++)
	{
		j=0;

		im =(short *)read_analyze_image(imageList[i], &nx, &ny, &nz, &dx, &dy, &dz, &type, 0);
		nv = nx*ny*nz;

		if(im!=NULL)
		{
			for(int v=0; v<nv; v++)
			{
				if(mask[v]>0) 
					X[i*p + j++] = im[v];
			}

			free(im);
		}
	}

	return(X);
}

short *readDataMatrixShort_nifti(char **imageList, int n, int p, short *mask)
{
   short *im;
   nifti_1_header hdr;
   //int nx, ny, nz, nv, type;
   int nx, ny, nz, nv;
   short *X;
   int j;

	X = (short *)calloc(n*p, sizeof(short));

	for(int i=0; i<n; i++)
	{
		j=0;

		//im =(short *)read_analyze_image(imageList[i], &nx, &ny, &nz, &dx, &dy, &dz, &type, 0);
        im = (short *)read_nifti_image(imageList[i], &hdr);
        nx = hdr.dim[1];
        ny = hdr.dim[2];
        nz = hdr.dim[3];

		nv = nx*ny*nz;

		if(im!=NULL)
		{
			for(int v=0; v<nv; v++)
			{
				if(mask[v]>0) 
					X[i*p + j++] = im[v];
			}

			free(im);
		}
	}

	return(X);
}

void ijk2xyz(float *T, int nx, int ny, int nz, float dx, float dy, float dz)
{
	T[0] =dx;	T[1] =0.0; 	T[2] =0.0; 	T[3] =-dx*(nx-1.0)/2.0;
	T[4] =0.0; 	T[5] =dy; 	T[6] =0.0; 	T[7] =-dy*(ny-1.0)/2.0;
	T[8] =0.0; 	T[9] =0.0; 	T[10]=dz; 	T[11]=-dz*(nz-1.0)/2.0;
	T[12]=0.0; 	T[13]=0.0; 	T[14]=0.0; 	T[15]=1.0;

	return;
}

void ijk2xyz(float *T, DIM dim)
{
	T[0] =dim.dx;	T[1] =0.0; 	T[2] =0.0; 	T[3] =-dim.dx*(dim.nx-1.0)/2.0;
	T[4] =0.0; 	T[5] =dim.dy; 	T[6] =0.0; 	T[7] =-dim.dy*(dim.ny-1.0)/2.0;
	T[8] =0.0; 	T[9] =0.0; 	T[10]=dim.dz; 	T[11]=-dim.dz*(dim.nz-1.0)/2.0;
	T[12]=0.0; 	T[13]=0.0; 	T[14]=0.0; 	T[15]=1.0;

	return;
}

void xyz2ijk(float *T, int nx, int ny, int nz, float dx, float dy, float dz)
{
	if(dx==0.0 || dy==0.0 || dz==0.0)
	{
		printf("\nWarning: zero voxel dimensions encountered!\n\n");
		return;
	}

	T[0] =1.0/dx;  	T[1] =0.0; 	T[2] =0.0; 	T[3] =(nx-1.0)/2.0;
	T[4] =0.0; 	T[5] =1.0/dy;  	T[6] =0.0; 	T[7] =(ny-1.0)/2.0;
	T[8] =0.0; 	T[9] =0.0; 	T[10]=1.0/dz; 	T[11]=(nz-1.0)/2.0;
	T[12]=0.0; 	T[13]=0.0; 	T[14]=0.0; 	T[15]=1.0;

	return;
}

void xyz2ijk(float *T, DIM dim)
{
	if(dim.dx==0.0 || dim.dy==0.0 || dim.dz==0.0)
	{
		printf("\nWarning: zero voxel dimensions encountered!\n\n");
		return;
	}

	T[0] =1.0/dim.dx;  	T[1] =0.0; 	T[2] =0.0; 	T[3] =(dim.nx-1.0)/2.0;
	T[4] =0.0; 	T[5] =1.0/dim.dy;  	T[6] =0.0; 	T[7] =(dim.ny-1.0)/2.0;
	T[8] =0.0; 	T[9] =0.0; 	T[10]=1.0/dim.dz; 	T[11]=(dim.nz-1.0)/2.0;
	T[12]=0.0; 	T[13]=0.0; 	T[14]=0.0; 	T[15]=1.0;

	return;
}

void orientationCodeConverter(int integerCode, char *stringCode)
{
	switch(integerCode)
	{
		case 1:
			stringCode[0]='P'; stringCode[1]='I';  stringCode[2]='L';
			break;
		case 2:
			stringCode[0]='P'; stringCode[1]='I';  stringCode[2]='R';
			break;
		case 3:
			stringCode[0]='P'; stringCode[1]='S';  stringCode[2]='L';
			break;
		case 4:
			stringCode[0]='P'; stringCode[1]='S';  stringCode[2]='R';
			break;
		case 5:
			stringCode[0]='P'; stringCode[1]='L';  stringCode[2]='I';
			break;
		case 6:
			stringCode[0]='P'; stringCode[1]='L';  stringCode[2]='S';
			break;
		case 7:
			stringCode[0]='P'; stringCode[1]='R';  stringCode[2]='I';
			break;
		case 8:
			stringCode[0]='P'; stringCode[1]='R';  stringCode[2]='S';
			break;

		case 9:
			stringCode[0]='A'; stringCode[1]='I';  stringCode[2]='L';
			break;
		case 10:
			stringCode[0]='A'; stringCode[1]='I';  stringCode[2]='R';
			break;
		case 11:
			stringCode[0]='A'; stringCode[1]='S';  stringCode[2]='L';
			break;
		case 12:
			stringCode[0]='A'; stringCode[1]='S';  stringCode[2]='R';
			break;
		case 13:
			stringCode[0]='A'; stringCode[1]='L';  stringCode[2]='I';
			break;
		case 14:
			stringCode[0]='A'; stringCode[1]='L';  stringCode[2]='S';
			break;
		case 15:
			stringCode[0]='A'; stringCode[1]='R';  stringCode[2]='I';
			break;
		case 16:
			stringCode[0]='A'; stringCode[1]='R';  stringCode[2]='S';
			break;

		case 17:
			stringCode[0]='I'; stringCode[1]='P';  stringCode[2]='L';
			break;
		case 18:
			stringCode[0]='I'; stringCode[1]='P';  stringCode[2]='R';
			break;
		case 19:
			stringCode[0]='I'; stringCode[1]='A';  stringCode[2]='L';
			break;
		case 20:
			stringCode[0]='I'; stringCode[1]='A';  stringCode[2]='R';
			break;
		case 21:
			stringCode[0]='I'; stringCode[1]='L';  stringCode[2]='P';
			break;
		case 22:
			stringCode[0]='I'; stringCode[1]='L';  stringCode[2]='A';
			break;
		case 23:
			stringCode[0]='I'; stringCode[1]='R';  stringCode[2]='P';
			break;
		case 24:
			stringCode[0]='I'; stringCode[1]='R';  stringCode[2]='A';
			break;

		case 25:
			stringCode[0]='S'; stringCode[1]='P';  stringCode[2]='L';
			break;
		case 26:
			stringCode[0]='S'; stringCode[1]='P';  stringCode[2]='R';
			break;
		case 27:
			stringCode[0]='S'; stringCode[1]='A';  stringCode[2]='L';
			break;
		case 28:
			stringCode[0]='S'; stringCode[1]='A';  stringCode[2]='R';
			break;
		case 29:
			stringCode[0]='S'; stringCode[1]='L';  stringCode[2]='P';
			break;
		case 30:
			stringCode[0]='S'; stringCode[1]='L';  stringCode[2]='A';
			break;
		case 31:
			stringCode[0]='S'; stringCode[1]='R';  stringCode[2]='P';
			break;
		case 32:
			stringCode[0]='S'; stringCode[1]='R';  stringCode[2]='A';
			break;

		case 33:
			stringCode[0]='L'; stringCode[1]='P';  stringCode[2]='I';
			break;
		case 34:
			stringCode[0]='L'; stringCode[1]='P';  stringCode[2]='S';
			break;
		case 35:
			stringCode[0]='L'; stringCode[1]='A';  stringCode[2]='I';
			break;
		case 36:
			stringCode[0]='L'; stringCode[1]='A';  stringCode[2]='S';
			break;
		case 37:
			stringCode[0]='L'; stringCode[1]='I';  stringCode[2]='P';
			break;
		case 38:
			stringCode[0]='L'; stringCode[1]='I';  stringCode[2]='A';
			break;
		case 39:
			stringCode[0]='L'; stringCode[1]='S';  stringCode[2]='P';
			break;
		case 40:
			stringCode[0]='L'; stringCode[1]='S';  stringCode[2]='A';
			break;

		case 41:
			stringCode[0]='R'; stringCode[1]='P';  stringCode[2]='I';
			break;
		case 42:
			stringCode[0]='R'; stringCode[1]='P';  stringCode[2]='S';
			break;
		case 43:
			stringCode[0]='R'; stringCode[1]='A';  stringCode[2]='I';
			break;
		case 44:
			stringCode[0]='R'; stringCode[1]='A';  stringCode[2]='S';
			break;
		case 45:
			stringCode[0]='R'; stringCode[1]='I';  stringCode[2]='P';
			break;
		case 46:
			stringCode[0]='R'; stringCode[1]='I';  stringCode[2]='A';
			break;
		case 47:
			stringCode[0]='R'; stringCode[1]='S';  stringCode[2]='P';
			break;
		case 48:
			stringCode[0]='R'; stringCode[1]='S';  stringCode[2]='A';
			break;
	}
}

int orientationCodeOK(char *stringCode)
{
	stringCode[0]=(char)toupper((int)(stringCode[0]));
	stringCode[1]=(char)toupper((int)(stringCode[1]));
	stringCode[2]=(char)toupper((int)(stringCode[2]));

	if(stringCode[0]=='P' && stringCode[1]=='I' &&  stringCode[2]=='L')
		return(1);
	if(stringCode[0]=='P' && stringCode[1]=='I' &&  stringCode[2]=='R')
		return(1);
	if(stringCode[0]=='P' && stringCode[1]=='S' &&  stringCode[2]=='L')
		return(1);
	if(stringCode[0]=='P' && stringCode[1]=='S' &&  stringCode[2]=='R')
		return(1);
	if(stringCode[0]=='P' && stringCode[1]=='L' &&  stringCode[2]=='I')
		return(1);
	if(stringCode[0]=='P' && stringCode[1]=='L' &&  stringCode[2]=='S')
		return(1);
	if(stringCode[0]=='P' && stringCode[1]=='R' &&  stringCode[2]=='I')
		return(1);
	if(stringCode[0]=='P' && stringCode[1]=='R' &&  stringCode[2]=='S')
		return(1);
	if(stringCode[0]=='A' && stringCode[1]=='I' &&  stringCode[2]=='L')
		return(1);
	if(stringCode[0]=='A' && stringCode[1]=='I' &&  stringCode[2]=='R')
		return(1);
	if(stringCode[0]=='A' && stringCode[1]=='S' &&  stringCode[2]=='L')
		return(1);
	if(stringCode[0]=='A' && stringCode[1]=='S' &&  stringCode[2]=='R')
		return(1);
	if(stringCode[0]=='A' && stringCode[1]=='L' &&  stringCode[2]=='I')
		return(1);
	if(stringCode[0]=='A' && stringCode[1]=='L' &&  stringCode[2]=='S')
		return(1);
	if(stringCode[0]=='A' && stringCode[1]=='R' &&  stringCode[2]=='I')
		return(1);
	if(stringCode[0]=='A' && stringCode[1]=='R' &&  stringCode[2]=='S')
		return(1);
	if(stringCode[0]=='I' && stringCode[1]=='P' &&  stringCode[2]=='L')
		return(1);
	if(stringCode[0]=='I' && stringCode[1]=='P' &&  stringCode[2]=='R')
		return(1);
	if(stringCode[0]=='I' && stringCode[1]=='A' &&  stringCode[2]=='L')
		return(1);
	if(stringCode[0]=='I' && stringCode[1]=='A' &&  stringCode[2]=='R')
		return(1);
	if(stringCode[0]=='I' && stringCode[1]=='L' &&  stringCode[2]=='P')
		return(1);
	if(stringCode[0]=='I' && stringCode[1]=='L' &&  stringCode[2]=='A')
		return(1);
	if(stringCode[0]=='I' && stringCode[1]=='R' &&  stringCode[2]=='P')
		return(1);
	if(stringCode[0]=='I' && stringCode[1]=='R' &&  stringCode[2]=='A')
		return(1);
	if(stringCode[0]=='S' && stringCode[1]=='P' &&  stringCode[2]=='L')
		return(1);
	if(stringCode[0]=='S' && stringCode[1]=='P' &&  stringCode[2]=='R')
		return(1);
	if(stringCode[0]=='S' && stringCode[1]=='A' &&  stringCode[2]=='L')
		return(1);
	if(stringCode[0]=='S' && stringCode[1]=='A' &&  stringCode[2]=='R')
		return(1);
	if(stringCode[0]=='S' && stringCode[1]=='L' &&  stringCode[2]=='P')
		return(1);
	if(stringCode[0]=='S' && stringCode[1]=='L' &&  stringCode[2]=='A')
		return(1);
	if(stringCode[0]=='S' && stringCode[1]=='R' &&  stringCode[2]=='P')
		return(1);
	if(stringCode[0]=='S' && stringCode[1]=='R' &&  stringCode[2]=='A')
		return(1);
	if(stringCode[0]=='L' && stringCode[1]=='P' &&  stringCode[2]=='I')
		return(1);
	if(stringCode[0]=='L' && stringCode[1]=='P' &&  stringCode[2]=='S')
		return(1);
	if(stringCode[0]=='L' && stringCode[1]=='A' &&  stringCode[2]=='I')
		return(1);
	if(stringCode[0]=='L' && stringCode[1]=='A' &&  stringCode[2]=='S')
		return(1);
	if(stringCode[0]=='L' && stringCode[1]=='I' &&  stringCode[2]=='P')
		return(1);
	if(stringCode[0]=='L' && stringCode[1]=='I' &&  stringCode[2]=='A')
		return(1);
	if(stringCode[0]=='L' && stringCode[1]=='S' &&  stringCode[2]=='P')
		return(1);
	if(stringCode[0]=='L' && stringCode[1]=='S' &&  stringCode[2]=='A')
		return(1);
	if(stringCode[0]=='R' && stringCode[1]=='P' &&  stringCode[2]=='I')
		return(1);
	if(stringCode[0]=='R' && stringCode[1]=='P' &&  stringCode[2]=='S')
		return(1);
	if(stringCode[0]=='R' && stringCode[1]=='A' &&  stringCode[2]=='I')
		return(1);
	if(stringCode[0]=='R' && stringCode[1]=='A' &&  stringCode[2]=='S')
		return(1);
	if(stringCode[0]=='R' && stringCode[1]=='I' &&  stringCode[2]=='P')
		return(1);
	if(stringCode[0]=='R' && stringCode[1]=='I' &&  stringCode[2]=='A')
		return(1);
	if(stringCode[0]=='R' && stringCode[1]=='S' &&  stringCode[2]=='P')
		return(1);
	if(stringCode[0]=='R' && stringCode[1]=='S' &&  stringCode[2]=='A')
		return(1);

	return(0);
}

/********************************************************************
This function removes space, tab and newline charactors from the
input charactor string.
********************************************************************/
void remove_space(char *inp)
{
   int i,j;
   int L;

   L=strlen(inp);

   j=0;
   for(i=0;i<L;i++)
   {
      if(inp[i]!=' ' && inp[i]!='\t' && inp[i]!='\n')
         inp[j++]=inp[i];
   }
   inp[j]='\0';
}

void sobel_edge_x(short *in, float *out, int nx, int ny)
{
   float sx;

   for(int i=1; i<nx-1; i++)
   for(int j=1; j<ny-1; j++)
   {
      sx = 
      in[(j-1)*nx + i+1] + 2.0*in[j*nx + i+1] + in[(j+1)*nx + i+1] 
     -in[(j-1)*nx + i-1] - 2.0*in[j*nx + i-1] - in[(j+1)*nx + i-1];

      out[j*nx + i] = sx;
   }
}

void sobel_edge_y(short *in, float *out, int nx, int ny)
{
   float sy;

   for(int i=1; i<nx-1; i++)
   for(int j=1; j<ny-1; j++)
   {
      sy = 
      in[(j+1)*nx + i-1] + 2.0*in[(j+1)*nx + i] + in[(j+1)*nx + i+1] 
     -in[(j-1)*nx + i-1] - 2.0*in[(j-1)*nx + i] - in[(j-1)*nx + i+1];

      out[j*nx + i] = sy;
   }
}

void sobel_edge(short *in, float *out, int nx, int ny)
{
   float sx, sy;

   for(int i=1; i<nx-1; i++)
   for(int j=1; j<ny-1; j++)
   {
      sx = 
      in[(j-1)*nx + i+1] + 2.0*in[j*nx + i+1] + in[(j+1)*nx + i+1] 
     -in[(j-1)*nx + i-1] - 2.0*in[j*nx + i-1] - in[(j+1)*nx + i-1];

      sy = 
      in[(j+1)*nx + i-1] + 2.0*in[(j+1)*nx + i] + in[(j+1)*nx + i+1] 
     -in[(j-1)*nx + i-1] - 2.0*in[(j-1)*nx + i] - in[(j-1)*nx + i+1];

      out[j*nx + i] = sqrtf(sx*sx + sy*sy);
   }
}

void sobel_edge(short *in, short *out, int nx, int ny)
{
   float sx, sy;

   for(int i=1; i<nx-1; i++)
   for(int j=1; j<ny-1; j++)
   {
      sx = 
      in[(j-1)*nx + i+1] + 2.0*in[j*nx + i+1] + in[(j+1)*nx + i+1] 
     -in[(j-1)*nx + i-1] - 2.0*in[j*nx + i-1] - in[(j+1)*nx + i-1];

      sy = 
      in[(j+1)*nx + i-1] + 2.0*in[(j+1)*nx + i] + in[(j+1)*nx + i+1] 
     -in[(j-1)*nx + i-1] - 2.0*in[(j-1)*nx + i] - in[(j-1)*nx + i+1];

      out[j*nx + i] = (short)(sqrtf(sx*sx + sy*sy) + 0.5);
   }
}

