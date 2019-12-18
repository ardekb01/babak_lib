/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \
*  3dwarper.c                                                *
*  Copyright 2004 by Babak A. Ardekani                       *
*  ALL RIGHTS RESERVED.  No part of this program may be      *
*  used, transferred, or modified by any means without       *
*  prior written permission.  Making copies of any part      *
*  of this program for any purpose is a violation of         *
*  copyright laws.                                           *
\ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <sys/resource.h>
#include <volume.h>
#include <spm_analyze.h>
#include <babak_lib.h>
#include <sph.h>
#include <smooth.h>
#include <minmax.h>
#include <interpolator.h>
#include <stats.h>

#define YES 1
#define NO 0
#define MAXITER 10
#define MAXR 15

#define XMATRIXSIZE 255
#define YMATRIXSIZE 255
#define ZMATRIXSIZE 189
#define VOXELSIZE 1.0

//////////////////////////////////////////////////////////////////////////////////////////////////

// multi-resolution image dimensions
DIM dim1, dim2, dim4, dim8;

int opt;

static struct option options[] =
{
   {"-version", 0, 'V'},
   {"-Version", 0, 'V'},
   {"-affine", 0, 'a'},
   {"-a", 0, 'a'},
   {"-I", 0, 'I'},
   {"-hr", 0, 'H'},

   {"-iter8", 1, '8'},
   {"-iter4", 1, '4'},
   {"-iter2", 1, '2'},
   {"-iter1", 1, '1'},

   {"-r", 1, 'r'},
   {"-R", 1, 'R'},

   {"-T", 1, 'T'},
   {"-sub", 1, 's'},
   {"-trg", 1, 't'},

   {"-trglm", 1, 'm'},
   {"-sublm", 1, 'M'},

   {"-v", 0, 'v'},
   {"-verbose", 0, 'v'},
   {"-h",0,'h'},
   {"-help",0,'h'},
   {"-u", 1, 'u'},
   {"-o", 1, 'o'},

   {"-cubicspline", 0, 'c'},
   {0, 0, 0}
};

int opt_cubicspline=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

extern float detect_lm(SPH &searchsph, SPH &testsph, SHORTIM testim, int lmcm[], SPH &refsph, int lm[]);

void print_help_and_exit()
{
   printf("\n\nUsage:\n"
   "\t3dwarper [-v or -verbose] [-h or -help] [-affine or -a] [-I] [-version]\n"
   "\t[-u <filename>] [-o <filename>] [-cubicspline] [-r patch radius] [-R search radius]\n"
   "\t-sub <subject image> -trg <target image>\n\n" 

   "Required arguments:\n\n"
   "\t-sub <subject image>\n"
   "\t\tSpecifies the image to be registred to the <target image>. <subject image> is expected\n"
   "\t\tto be of type short in NIFTI format.\n\n"

   "\t-trg <target image>\n"
   "\t\tSpecifies the image to which the <subject image> is to be registred. <target image> is\n"
   "\t\texpected to be of type short in NIFTI format.\n\n" 

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

   "\t-I\n"
   "\t\tDoes not perform initial linear registration of any kind.\n\n"

   "\t-affine or -a\n"
   "\t\tLimits the registration to an affine transformation.\n\n"

   "\t-u <filename>\n"
   "\t\tStores the displacement vector field in the specified <filename>\n"
   "\t\t(default: <filename>=<subject image>_wrp.nii).'\n\n"

   "\t-o <filename>\n"
   "\t\tStores the non-linearly transformed (registered) <subject image> in\n"
   "\t\tthe specified <filename> (default: <filename>=C<subject image>.nii).\n\n"

   "\t-cubicspline\n"
   "\t\tThe output transformed (registered) <subject image> is generated using the\n"
   "\t\tcubic spline interpolation (default is trilinear interpolation).\n\n"
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

void nlReg(short *sub, short *trg, short *msk, DIM dim, int r, int R, float *Xwarp, float *Ywarp, float *Zwarp, int S)
{
   float I[16];
   float *Xw, *Yw, *Zw;

   SPH subsph(r);
   SPH trgsph(r);
   SPH searchsph(R);

   SHORTIM subim;
   set_dim(subim,dim);
   subim.v=sub;

   SHORTIM trgim;
   set_dim(trgim,dim);
   trgim.v=trg;

   int v; // voxel index going from 0 to dim.nv-1
   int P[3], Q[3];

   Xw = (float *)calloc(dim.nv, sizeof(float));
   Yw = (float *)calloc(dim.nv, sizeof(float));
   Zw = (float *)calloc(dim.nv, sizeof(float));

   for(int i=0; i<dim.nv; i++) if(msk[i]==0) sub[i]=trg[i]=0;
    
   for(int k=0; k<dim.nz; k++)
   for(int j=0; j<dim.ny; j++)
   for(int i=0; i<dim.nx; i++)
   {
      v = k*dim.np + j*dim.nx + i;

      if(msk[v]<=0 || trg[v]<=0)
      {
         continue;
      }

      trgsph.set(trgim,i,j,k);
      standardize(trgsph.v,trgsph.v,trgsph.n);

      P[0]=i; P[1]=j; P[2]=k;
      detect_lm(searchsph, subsph, subim, P, trgsph, Q);

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

         if(trg[v]<=0)
         {
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

// finds an affine registration in "A" to takes points from sub to trg space
void affineReg(short *sub, short *trg, short *msk, DIM dim, int r, int R, float *A)
{
   SPH subsph(r);
   SPH trgsph(r);
   SPH searchsph(R);

   SHORTIM subim;
   set_dim(subim,dim);
   subim.v=sub;

   SHORTIM trgim;
   set_dim(trgim,dim);
   trgim.v=trg;

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
      v = k*dim.np+ j*dim.nx+ i;

      if(msk[v]<=0 || sub[v]<=0 || trg[v]<=0)
      {
         flg[v]=0;
         continue;
      }

      subsph.set(subim,i,j,k);
      standardize(subsph.v,subsph.n);

      PT[3*c + 0] = QT[3*c + 0] = i;
      PT[3*c + 1] = QT[3*c + 1] = j;
      PT[3*c + 2] = QT[3*c + 2] = k;
      detect_lm(searchsph, trgsph, trgim, PT+3*c, subsph, QT+3*c);

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

void generateMultiResolution(short *sub, DIM subdim, float *T, float *Xwarp, float *Ywarp, float *Zwarp, short *im1, short *im2, short *im4, short *im8)
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

   tmp=computeReslicedImage2(sub, subdim, dim1, Xw, Yw, Zw);
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
void combine_warps_and_trans(float *trgTPIL, float *Xwarp, float *Ywarp, float *Zwarp, float *Xout, float *Yout, float *Zout, DIM dim_trg)
{
   int v;
   float xc,yc,zc;
   float xc1,yc1,zc1;
   float x,y,z;   
   float x1,y1,z1;   
   float i1,j1,k1;   

   xc=dim_trg.dx*(dim_trg.nx-1.0)/2.0;     /* +---+---+ */
   yc=dim_trg.dy*(dim_trg.ny-1.0)/2.0;
   zc=dim_trg.dz*(dim_trg.nz-1.0)/2.0;

   xc1=dim1.dx*(dim1.nx-1.0)/2.0;     /* +---+---+ */
   yc1=dim1.dy*(dim1.ny-1.0)/2.0;
   zc1=dim1.dz*(dim1.nz-1.0)/2.0;
 
   for(int k=0;k<dim_trg.nz;k++) 
   for(int j=0;j<dim_trg.ny;j++) 
   for(int i=0;i<dim_trg.nx;i++) 
   {
      v = k*dim_trg.np + j*dim_trg.nx + i;

      // (i*dx-xc) converts from image coordinates (i,j,z) to (x,y,z) coordinates
      x = (i*dim_trg.dx - xc);
      y = (j*dim_trg.dy - yc);
      z = (k*dim_trg.dz - zc);

      x1 =  trgTPIL[0]*x +trgTPIL[1]*y +trgTPIL[2]*z  +trgTPIL[3];
      y1 =  trgTPIL[4]*x +trgTPIL[5]*y +trgTPIL[6]*z  +trgTPIL[7];
      z1 =  trgTPIL[8]*x +trgTPIL[9]*y +trgTPIL[10]*z +trgTPIL[11];

      i1 = (x1 + xc1) / dim1.dx;
      j1 = (y1 + yc1) / dim1.dy;
      k1 = (z1 + zc1) / dim1.dz;

      x1 += linearInterpolator(i1,j1,k1,Xwarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      y1 += linearInterpolator(i1,j1,k1,Ywarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      z1 += linearInterpolator(i1,j1,k1,Zwarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      
      Xout[v] = x1 - i*dim_trg.dx + xc;
      Yout[v] = y1 - j*dim_trg.dy + yc;
      Zout[v] = z1 - k*dim_trg.dz + zc;
   }
}

int main(int argc, char **argv)
{
   // input "subject image" full path
   char subjectImageFile[1024]=""; 
   // input "trget image" full path
   char targetImageFile[1024]="";

   // extracted subject image filename only without suffix
   char subprefix[1024]=""; 
   // extracted target image filename only without suffix
   char trgprefix[1024]="";

   // original subject and target images loaded 
   // from subjectImageFile and targetImageFile
   short *sub, *trg;

   // NIFTI headers for input subject and target images
   nifti_1_header sub_hdr, trg_hdr;

   // structures for holding subject and target image dimensions
   DIM dim_sub, dim_trg;

   // If this flag is set to YES, initial rigid-body + affine transformations are NOT performed
   // subTPIL, trgTPIL, and sub_to_trg are all set to the identity matrix
   char opt_I=NO;

   // a linear transformation that brings the subject/target image from its native space to PIL space
   float subTPIL[16], trgTPIL[16];

   // filenames where AC, PC, and VSPS are manually specific for subject/target image
   char sublmfile[1024]="";
   char trglmfile[1024]="";

   // number if iterations of NL registration at each resolution level
   int iter8=4;
   int iter4=3;
   int iter2=2;
   int iter1=1;

   // Subject and target PIL versions and brain mask at different resolutons
   short *subPIL1=NULL, *subPIL2=NULL, *subPIL4=NULL, *subPIL8=NULL;
   short *trgPIL1=NULL, *trgPIL2=NULL, *trgPIL4=NULL, *trgPIL8=NULL;
   short *mskPIL1=NULL, *mskPIL2=NULL, *mskPIL4=NULL, *mskPIL8=NULL;

   // 4x4 identity matrix
   float I[16];

   // space holder for a generic filename
   char filename[1024];

   int patch_r=3; // patch radius
   int search_R=3; // search radius 

   // subject to target image affine transformation
   float sub_to_trg[16];

   // if set to YES, only an affine registration is performed
   char opt_affine=NO;

   // affine transformation file specified at input
   char affineTransformationFile[1024]="";

   // if specified at the command line, the non-linearly transformed subject 
   // image will be saved under this name
   char outputfile[1024]="";

   // if specified at the command line, the final displacement field
   // image will be saved under this name
   char warpfile[1024]="";

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
         case 'I':
            opt_I= YES;
            break;
         case 'r':
            patch_r = atoi(optarg);
            break;
         case 'R':
            search_R = atoi(optarg);
            break;
         case 'a':
            opt_affine = YES;
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
         case 'T':
            sprintf(affineTransformationFile,"%s",optarg);
            break;
         case 'm':
            sprintf(trglmfile,"%s",optarg);
            break;
         case 'M':
            sprintf(sublmfile,"%s",optarg);
            break;
         case 'h':
            print_help_and_exit();
            break;
         case 'v':
            opt_v=YES;
            break;
         case 's':
            sprintf(subjectImageFile,"%s",optarg);
            break;
         case 't':
            sprintf(targetImageFile,"%s",optarg);
            break;
         case 'u':
            sprintf(warpfile,"%s",optarg);
            break;
         case 'o':
            sprintf(outputfile,"%s",optarg);
            break;
         case 'c':
            opt_cubicspline=YES;
            break;
         case '?':
            print_help_and_exit();
		}
	}

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Variable initializations
   getARTHOME();

   // Initialize dim1, dim2, dim4, and dim8 multi-resolution image dimensions
   dim1.nx = XMATRIXSIZE;
   dim1.ny = YMATRIXSIZE;
   dim1.nz = ZMATRIXSIZE;
   dim1.np=dim1.nx*dim1.ny; 
   dim1.nv=dim1.np*dim1.nz; 
   dim1.dx = dim1.dy = dim1.dz = VOXELSIZE;

   dim2.nx = XMATRIXSIZE/2;
   dim2.ny = YMATRIXSIZE/2;
   dim2.nz = ZMATRIXSIZE/2;
   dim2.np=dim2.nx*dim2.ny; 
   dim2.nv=dim2.np*dim2.nz; 
   dim2.dx = dim2.dy = dim2.dz = VOXELSIZE*2.0;

   dim4.nx = XMATRIXSIZE/4;
   dim4.ny = YMATRIXSIZE/4;
   dim4.nz = ZMATRIXSIZE/4;
   dim4.np=dim4.nx*dim4.ny; 
   dim4.nv=dim4.np*dim4.nz; 
   dim4.dx = dim4.dy = dim4.dz = VOXELSIZE*4.0;

   dim8.nx = XMATRIXSIZE/8;
   dim8.ny = YMATRIXSIZE/8;
   dim8.nz = ZMATRIXSIZE/8;
   dim8.np=dim8.nx*dim8.nx; 
   dim8.nv=dim8.np*dim8.nz; 
   dim8.dx = dim8.dy = dim8.dz = VOXELSIZE*8.0;

   trgPIL1 = (short *)calloc(dim1.nv,sizeof(short));
   subPIL1 = (short *)calloc(dim1.nv,sizeof(short));
   mskPIL1 = (short *)calloc(dim1.nv,sizeof(short));
   trgPIL2 = (short *)calloc(dim2.nv,sizeof(short));
   subPIL2 = (short *)calloc(dim2.nv,sizeof(short));
   mskPIL2 = (short *)calloc(dim2.nv,sizeof(short));
   trgPIL4 = (short *)calloc(dim4.nv,sizeof(short));
   subPIL4 = (short *)calloc(dim4.nv,sizeof(short));
   mskPIL4 = (short *)calloc(dim4.nv,sizeof(short));
   trgPIL8 = (short *)calloc(dim8.nv,sizeof(short));
   subPIL8 = (short *)calloc(dim8.nv,sizeof(short));
   mskPIL8 = (short *)calloc(dim8.nv,sizeof(short));

   // reset unreasonable values to default
   if(iter8 < 0 || iter8>MAXITER ) iter8=4;
   if(iter4 < 0 || iter4>MAXITER ) iter4=3;
   if(iter2 < 0 || iter2>MAXITER ) iter2=2;
   if(iter1 < 0 || iter1>MAXITER ) iter1=1;

   if(patch_r <= 0 || patch_r > MAXR ) patch_r=3;
   if(search_R <= 0 || search_R > MAXR ) search_R=3;

   ////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////
   // ensure that subject and target images are specified at the command line
   if( subjectImageFile[0] == '\0')
   {
      printf("\nPlease specify a \"subject image\" using -sub <filename.nii> ...\n\n");		
      exit(0);
   } 

   if( targetImageFile[0] == '\0')
   {
      printf("\nPlease specify a \"target image\" -trg <filename.nii> ...\n\n");		
      exit(0);
   }
   ////////////////////////////////////////////////////////////////////////////////////////////

   // extract subject and target filenames without path/suffix
   if( niftiFilename(subprefix, subjectImageFile)==0 ) { exit(0); }
   if( niftiFilename(trgprefix, targetImageFile)==0 ) { exit(0); }

   ////////////////////////////////////////////////////////////////////////////////////////////
   // read subject and target images
   sub = (int2 *)read_nifti_image(subjectImageFile, &sub_hdr);
   if(sub==NULL) 
   {
      printf("\nError: Reading subject image %s failed.\n\n",subjectImageFile);
      exit(0);
   }
   set_dim(dim_sub,sub_hdr); // transfer info from sub_hdr to dim_sub

   trg = (int2 *)read_nifti_image(targetImageFile, &trg_hdr);
   if(trg==NULL) 
   {
      printf("\nError: Reading target image %s failed.\n\n",targetImageFile);
      exit(0);
   }
   set_dim(dim_trg,trg_hdr); // transfer info from trg_hdr to dim_trg
   ////////////////////////////////////////////////////////////////////////////////////////////

   // print some info
   if(opt_v)
   {
      printf("Subject image file: %s\n",subjectImageFile);
      printf("\tMatrix size: %d x %d x %d\n",dim_sub.nx, dim_sub.ny, dim_sub.nz);
      printf("\tVoxel size: %6.4f x %6.4f x %6.4f\n",dim_sub.dx, dim_sub.dy, dim_sub.dz);
      printf("Target image file: %s\n",targetImageFile);
      printf("\tMatrix size: %d x %d x %d\n",dim_trg.nx, dim_trg.ny, dim_trg.nz);
      printf("\tVoxel size: %6.4f x %6.4f x %6.4f\n",dim_trg.dx, dim_trg.dy, dim_trg.dz);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////
   // find subTPIL and trgTPIL using automated MSP, AC/PC and 8 MSP landmarks detection
   if(opt_I==YES)
   {
      set_to_I(subTPIL,4);
      set_to_I(trgTPIL,4);
   }
   else {
      // find subTPIL 
      char orient[4]="";

      if(opt_v) printf("Computing subject image PIL transformation ...\n");
      if(sublmfile[0] != '\0' && opt_v) printf("Subject image landmarks are read from %s\n",sublmfile);
      orient[0]='\0';
      new_PIL_transform(subjectImageFile, sublmfile, orient, subTPIL, 1);

      // find trgTPIL
      if(opt_v) printf("Computing target image PIL transformation ...\n");
      if(trglmfile[0] != '\0' && opt_v) printf("Target image landmarks are read from %s\n",trglmfile);
      orient[0]='\0';
      new_PIL_transform(targetImageFile, trglmfile, orient, trgTPIL, 1);
   }
   ////////////////////////////////////////////////////////////////////////////////////////////
   
   // compute trgPIL1, trgPIL2, trgPIL4 and trgPIL8
   generateMultiResolution(trg, dim_trg, trgTPIL, trgPIL1, trgPIL2, trgPIL4, trgPIL8);

   ////////////////////////////////////////////////////////////////////////////////////////////
   // compute mskPIL1, mskPIL2, mskPIL4 and mskPIL8
   if(opt_I == YES)
   {
      set_to_I(I,4);
      generateMultiResolution(sub, dim_sub, I, mskPIL1, mskPIL2, mskPIL4, mskPIL8);
   }
   else
   {
      int2 *PILbraincloud;
      DIM PILbraincloud_dim;
      nifti_1_header PILbraincloud_hdr;

      sprintf(filename,"%s/PILbrain.nii",ARTHOME);

      PILbraincloud = (int2 *)read_nifti_image(filename, &PILbraincloud_hdr);

      if(PILbraincloud==NULL)
      {
            printf("Error reading %s, aborting ...\n", filename);
            exit(1);
      }

      set_dim(PILbraincloud_dim, PILbraincloud_hdr);

      set_to_I(I,4);
      generateMultiResolution(PILbraincloud, PILbraincloud_dim, I, mskPIL1, mskPIL2, mskPIL4, mskPIL8);

      free(PILbraincloud);
   }
   ////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Compute the affine transformation  sub_to_trg and update subTPIL

   if(opt_v) printf("Patch radius = %d mm\n", patch_r);
   if(opt_v) printf("Search radius = %d mm\n", search_R);

   if(opt_I == YES)
   {
      set_to_I(sub_to_trg,4);
   }
   else if( affineTransformationFile[0] != '\0' ) // -I option trumps -T
   {
      if( loadTransformation(affineTransformationFile,sub_to_trg)== 0 ) // no error
         multi(trgTPIL,4,4,sub_to_trg,4,4,subTPIL);  // update subTPIL according to A
      else
      {
         printf("Reading %s failed, aborting ...\n",affineTransformationFile);
         exit(0);
      }
   }
   else
   {
      float *invT;		
      float A[16];
      FILE *fp;

      if(opt_v) printf("Affine registration @ 12.5%% resolution ...\n");
      generateMultiResolution(sub, dim_sub, subTPIL, subPIL1, subPIL2, subPIL4, subPIL8);
      affineReg(subPIL8, trgPIL8, mskPIL8, dim8, patch_r, search_R, A);
      multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

      if(opt_v) printf("Affine registration @ 25%% resolution ...\n");
      generateMultiResolution(sub, dim_sub, subTPIL, subPIL1, subPIL2, subPIL4, subPIL8);
      affineReg(subPIL4, trgPIL4, mskPIL4, dim4, patch_r, search_R, A);
      multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

      if(opt_v) printf("Affine registration @ 50%% resolution ...\n");
      generateMultiResolution(sub, dim_sub, subTPIL, subPIL1, subPIL2, subPIL4, subPIL8);
      affineReg(subPIL2, trgPIL2, mskPIL2, dim2, patch_r, search_R, A);
      multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

      invT = inv4(trgTPIL);
      multi(invT, 4,4, subTPIL, 4, 4, sub_to_trg);
      free(invT);

      sprintf(filename,"%s_affine.mrx",subprefix);
      if(opt_v) printf("Saving subject image -> target image affine transformaion in %s ...\n",filename);
      fp = fopen(filename,"w");
      if(fp==NULL) file_open_error(filename);
      printMatrix(sub_to_trg,4,4,"",fp);
      fclose(fp);
   }

   if(opt_v) printMatrix(sub_to_trg,4,4,"subject image -> target image affine transformation",NULL);
   ////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Save affine transformed subject image
   {
      short *tmp;
      float *invT;		

      invT=inv4(sub_to_trg);
      tmp = resliceImage(sub, dim_sub, dim_trg, invT, LIN); 
      free(invT);

      sprintf(filename,"%s_affine.nii",subprefix);
      if(opt_v) printf("Saving subject image affine transformed to target space in %s ...\n",filename);
      save_nifti_image(filename, tmp, &trg_hdr);
      free(tmp);
   }

   if(opt_affine) 
   {
      if(opt_v) printf("THE END\n");
      exit(0);
   }
   //////////////////////////////////////////////////////////////////
   
   //////////////////////////////////////////////////////////////////
   // Non-linear registration from trgPIL to subPIL
   
   // intermediate displacement field from trgPIL to subPIL
   float *Xwarp=NULL, *Ywarp=NULL, *Zwarp=NULL;

   Xwarp = (float *)calloc(dim1.nv, sizeof(float));
   Ywarp = (float *)calloc(dim1.nv, sizeof(float));
   Zwarp = (float *)calloc(dim1.nv, sizeof(float));

   if(iter8>=1)
   {
      if(opt_v) printf("Non-linear registration @ 12.5%% resolution ...\n");
      for(int i=0; i<iter8; i++)
      {
         generateMultiResolution(sub, dim_sub, subTPIL, Xwarp, Ywarp, Zwarp, subPIL1, subPIL2, subPIL4, subPIL8);
         nlReg(subPIL8, trgPIL8, mskPIL8, dim8, patch_r, search_R/3, Xwarp, Ywarp, Zwarp,5);
      }
   }

   if(iter4>=1)
   {
      if(opt_v) printf("Non-linear registration @ 25%% resolution ...\n");
      for(int i=0; i<iter4; i++)
      {
         generateMultiResolution(sub, dim_sub, subTPIL, Xwarp, Ywarp, Zwarp, subPIL1, subPIL2, subPIL4, subPIL8);
         nlReg(subPIL4, trgPIL4, mskPIL4, dim4, patch_r, 2*search_R/3, Xwarp, Ywarp, Zwarp,5);
      }
   }

   if(iter2>=1)
   {
      if(opt_v) printf("Non-linear registration @ 50%% resolution ...\n");
      for(int i=0; i<iter2; i++)
      {
         generateMultiResolution(sub, dim_sub, subTPIL, Xwarp, Ywarp, Zwarp, subPIL1, subPIL2, subPIL4, subPIL8);
         nlReg(subPIL2, trgPIL2, mskPIL2, dim2, patch_r, search_R, Xwarp, Ywarp, Zwarp,5);
      }
   }

   if(iter1>=1)
   {
      if(opt_v) printf("Non-linear registration @ 100%% resolution ...\n");
      for(int i=0; i<iter1; i++)
      {
         generateMultiResolution(sub, dim_sub, subTPIL, Xwarp, Ywarp, Zwarp, subPIL1, subPIL2, subPIL4, subPIL8);
         nlReg(subPIL1, trgPIL1, mskPIL1, dim1, patch_r, search_R, Xwarp, Ywarp, Zwarp,5);
      }
   }
   //////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////
   // combine trgTPIL, (Xwarp,Ywarp,Zwarp) and subTPIL to obtain (Xout,Yout,Zout)

   // final displacement field from trg to sub
   float *Xout=NULL, *Yout=NULL, *Zout=NULL; 

   Xout = (float *)calloc(dim_trg.nv, sizeof(float));
   Yout = (float *)calloc(dim_trg.nv, sizeof(float));
   Zout = (float *)calloc(dim_trg.nv, sizeof(float));

   combine_warps_and_trans(dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim1.dy, dim1.dz, Xwarp, Ywarp, Zwarp, subTPIL);
   combine_warps_and_trans(trgTPIL, Xwarp, Ywarp, Zwarp, Xout, Yout, Zout, dim_trg);

   ////////////////////////////////////////////////////////////////////////////////////////////
   // saved warped subject image
   {
      short *Csub;

      Csub=computeReslicedImage2(sub, dim_sub, dim_trg, Xout, Yout, Zout);

      if(outputfile[0]=='\0') sprintf(filename,"C%s.nii",subprefix);
      else sprintf(filename,"%s",outputfile);
      if(opt_v) printf("Saving non-linearly transformed subject image in %s ...\n",filename);
      save_nifti_image(filename, Csub, &trg_hdr);

      free(Csub);
   }
   ////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////////////////
   // save (Xout, Yout, Zout)
   {
      float min=0.0, max=0.0, s=0.0;
      short *sdum;
      nifti1_extender extender;
      FILE *fp;

      if(warpfile[0]=='\0') sprintf(warpfile,"%s_wrp.nii",subprefix);

      if(opt_v) printf("Saving displacement field %s ...\n",warpfile);

      minmax(Xout, dim_trg.nv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Yout, dim_trg.nv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Zout, dim_trg.nv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      fp=fopen(warpfile,"w");
      if(fp==NULL) file_open_error(warpfile);
      trg_hdr.bitpix=8*sizeof(short);
      trg_hdr.vox_offset=352.0;
      trg_hdr.dim[0]=5;
      trg_hdr.dim[4]=1;
      trg_hdr.dim[5]=3;
      trg_hdr.pixdim[4] = 0.0;
      trg_hdr.datatype=DT_SIGNED_SHORT;
      trg_hdr.intent_code=NIFTI_INTENT_VECTOR;
      trg_hdr.scl_slope=s/32000;
      extender.extension[0]=0;
      fwrite(&trg_hdr,sizeof(nifti_1_header),1,fp);
      fwrite(&extender, sizeof(nifti1_extender),1,fp);

      sdum = (short *)calloc(dim_trg.nv, sizeof(short));
      for(int i=0; i<dim_trg.nv; i++)
      {
         sdum[i] = (short)(Xout[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),dim_trg.nv,fp);

      for(int i=0; i<dim_trg.nv; i++)
      {
         sdum[i] = (short)(Yout[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),dim_trg.nv,fp);

      for(int i=0; i<dim_trg.nv; i++)
      {
         sdum[i] = (short)(Zout[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),dim_trg.nv,fp);

      fclose(fp);

      free(sdum);
   }
   ////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////////////////
   // free allocated memory 
   free(Xout); free(Yout); free(Zout);
   free(Xwarp); free(Ywarp); free(Zwarp);
   free(sub);
   free(trg);
   free(trgPIL1); free(trgPIL2); free(trgPIL4); free(trgPIL8);
   free(subPIL1); free(subPIL2); free(subPIL4); free(subPIL8);
   free(mskPIL1); free(mskPIL2); free(mskPIL4); free(mskPIL8);
   ////////////////////////////////////////////////////////////////////////////////////////////

   if(opt_v) printf("THE END\n");
}
