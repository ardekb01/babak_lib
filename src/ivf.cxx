#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <unistd.h>
#include <volume.h>
#include <spm_analyze.h>
#include <babak_lib.h>
#include <minmax.h>

#define YES 1
#define NO 0

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-v", 0, 'v'},
   {"-verbose", 0, 'v'},
   {"-h", 0, 'h'},
   {"-help", 0, 'h'},

   {"-FWHM", 1, 'F'},
   {"-nx", 1, 'x'},
   {"-ny", 1, 'y'},
   {"-nz", 1, 'z'},
   {"-dx", 1, 'X'},
   {"-dy", 1, 'Y'},
   {"-dz", 1, 'Z'},
   {"-i", 1, 'i'},
   {"-o", 1, 'o'},
   {0, 0, 0}
};

int opt_nx=NO;
int opt_ny=NO;
int opt_nz=NO;
int opt_dx=NO;
int opt_dy=NO;
int opt_dz=NO;
int opt_i=NO;
int opt_o=NO;

//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit()
{
   printf("\n\nUsage:\n"
   "\tivf [-v or -verbose] [-h or -help] [-FWHM <mm>] [-nx <columns> -ny <rows> -nz <slices>]\n"
   "\t[-dx <mm> -dy <mm> -dz <mm>] -i <input field> -o <output field>\n\n"

   "Required arguments:\n\n"
   "\t-i <input field>\n"
   "\t\tInput vector field in NIFTI format.\n\n"

   "\t-o <output field>\n"
   "\t\tName of the output vector field.\n\n"

   "Optional arguments:\n\n"
   "\t-v or -verbose\n"
   "\t\tEnables verbose mode.\n\n"

   "\t-h or -help\n"
   "\t\tPrints help message.\n\n"

   "\t-FWHM <mm>\n"
   "\t\tFull width at half maximum (FWHM) of the Gaussian kernel used in the inversion algorithm\n"
   "\t\t(default = 1.0 mm).\n\n"

   "\t-nx <columns>\n"
   "\t\tOutput vector field number of columns (default: same as input).\n\n"

   "\t-ny <rows>\n"
   "\t\tOutput vector field number of rows (default: same as input).\n\n"

   "\t-nz <slices>\n"
   "\t\tOutput vector field number of slices (default: same as input).\n\n"

   "\t-dx <mm>\n"
   "\t\tOutput vector field voxel dimension along x (increasing column-number) direction\n"
   "\t\t(default: same as input).\n\n"

   "\t-dy <mm>\n"
   "\t\tOutput vector field voxel dimension along y (increasing row-number) direction\n"
   "\t\t(default: same as input).\n\n"

   "\t-dz <mm>\n"
   "\t\tOutput vector field voxel dimension along z (increasing slice-number) direction\n"
   "\t\t(default: same as input).\n\n"

   );

   exit(0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

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
//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   int ii, jj, kk;
   float xx, yy, zz;

   int v1, v1_part1, v1_part2;
   int v2, v2_part1, v2_part2;

   float FWHM=1.0; 
   float sigma, K;
   int wx, wy, wz;
   float x,y,z;

   short *sdum;

   FILE *fp;

   nifti_1_header hdr;
   nifti1_extender extender;

   char filename2[1024];
   int nx2=0, ny2=0, nz2=0, np2=0, nv2=0;
   float dx2=0.0, dy2=0.0, dz2=0.0;
   float *Xwarp2, *Ywarp2, *Zwarp2;
   float x2, y2, z2;
   float xc2, yc2, zc2;

   char filename1[1024];
   int nx1, ny1, nz1, np1, nv1;
   float dx1, dy1, dz1;
   float *Xwarp1, *Ywarp1, *Zwarp1;
   float x1, y1, z1;
   int i2, j2, k2;
   float xc1, yc1, zc1;

   float *W,w;

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) 
      {
         case 'h':
            print_help_and_exit();
         case 'v':
            opt_v = YES;
            break;
         case 'F':
            FWHM=atof(optarg);
            break;
         case 'x':
            nx2=atoi(optarg);
            opt_nx = YES;
            break;
         case 'y':
            ny2=atoi(optarg);
            opt_ny = YES;
            break;
         case 'z':
            nz2=atoi(optarg);
            opt_nz = YES;
            break;
         case 'X':
            dx2=atof(optarg);
            opt_dx = YES;
            break;
         case 'Y':
            dy2=atof(optarg);
            opt_dy = YES;
            break;
         case 'Z':
            dz2=atof(optarg);
            opt_dz = YES;
            break;
         case 'i':
            sprintf(filename1,"%s",optarg);
            opt_i = YES;
            break;
         case 'o':
            sprintf(filename2,"%s",optarg);
            opt_o = YES;
            break;
         case '?':
            print_help_and_exit();
      }
   }

   ////////////////////////////////////////////////////////////////////////////////

   if(!opt_i)
   {
      printf("\nPlease specify the input vector field using the -i argument.\n\n");
      exit(0);
   }

   if(!opt_o)
   {
      printf("\nPlease specify a name for the output vector field using the -o argument.\n\n");
      exit(0);
   }
   ////////////////////////////////////////////////////////////////////////////////

   fp=fopen(filename1,"r");
   if(fp==NULL) file_open_error(filename1);

   fread(&hdr,sizeof(nifti_1_header),1,fp);
   fread(&extender, sizeof(nifti1_extender), 1, fp);

   nx1 = hdr.dim[1];
   ny1 = hdr.dim[2];
   nz1 = hdr.dim[3];
   dx1 = hdr.pixdim[1];
   dy1 = hdr.pixdim[2];
   dz1 = hdr.pixdim[3];
   nv1 = nx1 * ny1 * nz1;
   np1 = nx1 * ny1;

   sdum = (short *)calloc(nv1, sizeof(short));

   Xwarp1 = (float *)calloc(nv1,sizeof(float));
   Ywarp1 = (float *)calloc(nv1,sizeof(float));
   Zwarp1 = (float *)calloc(nv1,sizeof(float));

   fread(sdum,sizeof(short),nv1,fp);
   for(int i=0; i<nv1; i++) Xwarp1[i] = sdum[i]*hdr.scl_slope;

   fread(sdum,sizeof(short),nv1,fp);
   for(int i=0; i<nv1; i++) Ywarp1[i] = sdum[i]*hdr.scl_slope;

   fread(sdum,sizeof(short),nv1,fp);
   for(int i=0; i<nv1; i++) Zwarp1[i] = sdum[i]*hdr.scl_slope;

   delete sdum;

   fclose(fp);

   ////////////////////////////////////////////////////////////////////////////////

   if(!opt_nx || nx2<0) nx2=nx1;
   if(!opt_ny || ny2<0) ny2=ny1;
   if(!opt_nz || nz2<0) nz2=nz1;

   if(!opt_dx || dx2<0.0) dx2=dx1;
   if(!opt_dy || dy2<0.0) dy2=dy1;
   if(!opt_dz || dz2<0.0) dz2=dz1;

   ////////////////////////////////////////////////////////////////////////////////

   sigma = FWHM/2.35482; // converts FWHM to standard deviation
   K = 2.0*sigma*sigma;

   // Twice the standard deviation is greater than 95% of the area
   wx = (int)ceilf(2.0*sigma/dx2);
   wy = (int)ceilf(2.0*sigma/dy2);
   wz = (int)ceilf(2.0*sigma/dz2);

   if(opt_v)
   {
      printf("Standard deviation of the Gaussian kernel = %f\n\n",sigma);
      printf("Support region size: %d x %d x %d voxels\n\n", 2*wx+1, 2*wy+1, 2*wz+1);
   }

   ////////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////////////

   if(opt_v)
   {
      printf("Input displacement field = %s\n", filename1);
      printf("Input matrix size = %d x %d x %d\n",nx1, ny1, nz1);
      printf("Input voxel size = %7.5f x %7.5f x %7.5f mm^3\n\n",dx1, dy1, dz1);

      printf("Output displacement field = %s\n", filename2);
      printf("Output matrix size = %d x %d x %d\n",nx2, ny2, nz2);
      printf("Output voxel size = %7.5f x %7.5f x %7.5f mm^3\n",dx2, dy2, dz2);
   }

   ////////////////////////////////////////////////////////////////////////////////

   nv2 = nx2 * ny2 * nz2;
   np2 = nx2 * ny2;
   Xwarp2 = (float *)calloc(nv2,sizeof(float));
   Ywarp2 = (float *)calloc(nv2,sizeof(float));
   Zwarp2 = (float *)calloc(nv2,sizeof(float));

   W = (float *)calloc(nv2,sizeof(float));

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

   ////////////////////////////////////////////////////////////////////////////////
   {
      float min=0.0, max=0.0, s=0.0;

      minmax(Xwarp2, nv2, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Ywarp2, nv2, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Zwarp2, nv2, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      fp=fopen(filename2,"w");
      if(fp==NULL) file_open_error(filename2);
      hdr.dim[1]=nx2;
      hdr.dim[2]=ny2;
      hdr.dim[3]=nz2;
      hdr.pixdim[1]=dx2;
      hdr.pixdim[2]=dy2;
      hdr.pixdim[3]=dz2;
      hdr.bitpix=8*sizeof(short);
      hdr.vox_offset=352.0;
      hdr.dim[0]=5;
      hdr.dim[4]=1;
      hdr.dim[5]=3;
      hdr.pixdim[4] = 0.0;
      hdr.datatype=DT_SIGNED_SHORT;
      hdr.intent_code=NIFTI_INTENT_VECTOR;
      hdr.scl_slope=s/32000;
      extender.extension[0]=0;
      fwrite(&hdr,sizeof(nifti_1_header),1,fp);
      fwrite(&extender, sizeof(nifti1_extender),1,fp);

      sdum = (short *)calloc(nv2, sizeof(short));

      for(int i=0; i<nv2; i++)
      {
         sdum[i] = (short)(Xwarp2[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),nv2,fp);

      for(int i=0; i<nv2; i++)
      {
         sdum[i] = (short)(Ywarp2[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),nv2,fp);

      for(int i=0; i<nv2; i++)
      {
         sdum[i] = (short)(Zwarp2[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),nv2,fp);

      fclose(fp);

      free(sdum);
   }
}
