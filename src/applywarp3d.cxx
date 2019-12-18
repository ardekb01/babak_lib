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
#include <interpolator.h>

#define YES 1
#define NO 0

#define TYPEUNSIGNEDCHAR 2

short *computeReslicedImage(short *im1, nifti_1_header hdr1, nifti_1_header hdr2, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
float *computeReslicedImage(float *im1, nifti_1_header hdr1, nifti_1_header hdr2, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
unsigned char *computeReslicedImage(unsigned char *im1, nifti_1_header hdr1, nifti_1_header hdr2, float *Xwarp, float *Ywarp, float *Zwarp, float *T);

float *computeReslicedImage(float *im1, nifti_1_header hdr1, nifti_1_header hdr2, nifti_1_header hdr3,
float *Xwarp, float *Ywarp, float *Zwarp, float *T, float *Xw2d, float *Yw2d, float *Zw2d);

unsigned char *computeReslicedImage(unsigned char *im1, nifti_1_header hdr1, nifti_1_header hdr2, nifti_1_header hdr3,
float *Xwarp, float *Ywarp, float *Zwarp, float *T, float *Xw2d, float *Yw2d, float *Zw2d);

short *computeReslicedImage(short *im1, nifti_1_header hdr1, nifti_1_header hdr2, nifti_1_header hdr3,
float *Xwarp, float *Ywarp, float *Zwarp, float *T, float *Xw2d, float *Yw2d, float *Zw2d);

extern float *resizeXYZ(float *image1,
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-o", 1, 'o'},
   {"-nx", 1, 'x'},
   {"-ny", 1, 'y'},
   {"-nz", 1, 'z'},
   {"-dx", 1, 'X'},
   {"-dy", 1, 'Y'},
   {"-dz", 1, 'Z'},
        
   {"-h", 0, 'h'},
   {"-help", 0, 'h'},
   {"-version", 0, 'V'},
   {"-V", 0, 'V'},

   {"-secret", 0, 's'},
   {"-nn", 0, 'n'},
   {"-w", 1, '1'},
   {"-T", 1, '2'},
   {"-w2d", 1, '3'},
   {"-cubicspline", 0, '4'},
   {0, 0, 0}
};

int opt_o=NO;
int opt_nn=NO;
int opt_secret=NO;
int opt_w=NO;
int opt_T=NO;
int opt_w2d=NO;
int opt_cubicspline=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit();

int main(int argc, char **argv)
{
   int number_of_elements_read;

   nifti1_extender extender;
   nifti_1_header sub_hdr;
   nifti_1_header trg_hdr;
   int N;

	FILE *fp;

	// original target and object files stored by the 3dwarper program
   char prefix[1024]; // output file prefix

	int Tnx=0,Tny=0,Tnz=0,Tnv=0;

	float Tdx=0.0,Tdy=0.0,Tdz=0.0;

	float *Xwarp=NULL;
	float *Ywarp=NULL;
	float *Zwarp=NULL;

	float T[16];
	float *invT;	// inverse of matrix T

	char *dum;
	char outputfile[1024];
	char objfile[1024];
	char warpfile[1024];
	char Tfile[1024];

	char w2dfile[1024];

	char *op_im=NULL;

	dum = (char *)malloc(1024);

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'h':
            print_help_and_exit();
            exit(0);
         case 'V':
            printf("\nLast revision: July 12, 2012.\n\n");
            exit(0);
         case 'o':
            opt_o=YES;
            sprintf(prefix,"%s",optarg);
            break;
			case 'n':
				opt_nn=YES;
				break;
			case 's':
				opt_secret=YES;
				break;
			case '1':
				sprintf(warpfile,"%s",optarg);
				opt_w=YES;
				break;
			case '2':
				sprintf(Tfile,"%s",optarg);
				opt_T=YES;
				break;
			case '3':
				sprintf(w2dfile,"%s",optarg);
				opt_w2d=YES;
				break;
			case '4':
				opt_cubicspline=YES;
				break;
			case 'x':
				Tnx=atoi(optarg);
				break;
			case 'y':
				Tny=atoi(optarg);
				break;
			case 'z':
				Tnz=atoi(optarg);
				break;
			case 'X':
				Tdx=atof(optarg);
				break;
			case 'Y':
				Tdy=atof(optarg);
				break;
			case 'Z':
				Tdz=atof(optarg);
				break;
			case '?':
				print_help_and_exit();
		}
	}

////////////////////////////////////////////////////////////////////////////////////////////

   if( opt_o ) printf("\nOutput file prefix = %s\n",prefix);

////////////////////////////////////////////////////////////////////////////////////////////

   if( opt_T ) 
   {
      loadTransformation(Tfile, T);
   } 
   else
   {
      for(int i=0; i<16; i++) T[i]=0.0;
      T[0]=T[5]=T[10]=T[15]=1.0;
   }

   printf("\nInitial Transformation:");
   printf("\n%f\t%f\t%f\t%f",T[0],T[1],T[2],T[3]);
   printf("\n%f\t%f\t%f\t%f",T[4],T[5],T[6],T[7]);
   printf("\n%f\t%f\t%f\t%f",T[8],T[9],T[10],T[11]);
   printf("\n%f\t%f\t%f\t%f",T[12],T[13],T[14],T[15]);
   printf("\n");

//////////////////////////////////////////////////////////////////////////////////////////

   if(opt_w)
   {
      short *sdum;

      printf("\nDisplacement field = %s\n",warpfile);

      fp=fopen(warpfile,"r");
      if(fp==NULL)
      {
         printf("\nError opening %s, aborting ...\n", warpfile);
         exit(0);
      }

      number_of_elements_read = fread(&trg_hdr,sizeof(nifti_1_header),1,fp);
      if( number_of_elements_read != 1)
      {
         printf("\nError reading %s, aborting ...\n", warpfile);
         fclose(fp);
         exit(0);
      }

      number_of_elements_read = fread(&extender, sizeof(nifti1_extender), 1, fp);
      if( number_of_elements_read != 1)
      {
         printf("\nError reading %s, aborting ...\n", warpfile);
         fclose(fp);
         exit(0);
      }

      if(trg_hdr.dim[0]<1 || trg_hdr.dim[0]>7)
      {
         swapniftiheader(&trg_hdr);
      }
      Tnx = trg_hdr.dim[1];
      Tny = trg_hdr.dim[2];
      Tnz = trg_hdr.dim[3];
      Tdx = trg_hdr.pixdim[1];
      Tdy = trg_hdr.pixdim[2];
      Tdz = trg_hdr.pixdim[3];

      Tnv=Tnx*Tny*Tnz;

      printf("Target voxel size: %f %f %f\n",Tdx,Tdy,Tdz);
      printf("Target matrix size: %d %d %d\n",Tnx,Tny,Tnz);
      
      sdum = (short *)calloc(Tnv, sizeof(short));
      if( sdum == NULL )
      {
         fclose(fp);
         memory_allocation_error("sdum");
      }

      Xwarp=(float *)calloc(Tnv,sizeof(float));
      if( Xwarp == NULL )
      {
         fclose(fp);
         memory_allocation_error("Xwarp");
      }

      Ywarp=(float *)calloc(Tnv,sizeof(float));
      if( Ywarp == NULL )
      {
         fclose(fp);
         memory_allocation_error("Ywarp");
      }

      Zwarp=(float *)calloc(Tnv,sizeof(float));
      if( Zwarp == NULL )
      {
         fclose(fp);
         memory_allocation_error("Zwarp");
      }

      number_of_elements_read = fread(sdum,sizeof(short),Tnv,fp);
      if( number_of_elements_read != Tnv)
      {
         printf("\nError reading %s, aborting ...\n", warpfile);
         fclose(fp);
         exit(0);
      }
      for(int i=0; i<Tnv; i++) Xwarp[i] = sdum[i]*trg_hdr.scl_slope;

      number_of_elements_read = fread(sdum,sizeof(short),Tnv,fp);
      if( number_of_elements_read != Tnv)
      {
         printf("\nError reading %s, aborting ...\n", warpfile);
         fclose(fp);
         exit(0);
      }
      for(int i=0; i<Tnv; i++) Ywarp[i] = sdum[i]*trg_hdr.scl_slope;

      number_of_elements_read = fread(sdum,sizeof(short),Tnv,fp);
      if( number_of_elements_read != Tnv)
      {
         printf("\nError reading %s, aborting ...\n", warpfile);
         fclose(fp);
         exit(0);
      }
      for(int i=0; i<Tnv; i++) Zwarp[i] = sdum[i]*trg_hdr.scl_slope;

      fclose(fp);

      delete sdum;
   }
   else
   {
      trg_hdr.pixdim[1]=Tdx; trg_hdr.pixdim[2]=Tdy; trg_hdr.pixdim[3]=Tdz;
      trg_hdr.dim[1]=Tnx; trg_hdr.dim[2]=Tny; trg_hdr.dim[3]=Tnz;

      Tnv=Tnx*Tny*Tnz;
      Xwarp=(float *)calloc(Tnv,sizeof(float));
      if( Xwarp == NULL )
      {
         memory_allocation_error("Xwarp");
      }

      Ywarp=(float *)calloc(Tnv,sizeof(float));
      if( Ywarp == NULL )
      {
         memory_allocation_error("Ywarp");
      }

      Zwarp=(float *)calloc(Tnv,sizeof(float));
      if( Zwarp == NULL )
      {
         memory_allocation_error("Zwarp");
      }
   }

//////////////////////////////////////////////////////////////////////////////////////////

   N=argc-optind;

   invT=inv4(T);

   for(int i=0; i<N; i++)
   {
      char *sub_im;

      sprintf(objfile,"%s",argv[i+optind]);

      printf("\nTransforming image %s ...\n",objfile);

      if(opt_o)
      {
         if(N>1) sprintf(outputfile,"%s%03d.nii",prefix,i+1);
         else sprintf(outputfile,"%s.nii",prefix);
      }
      else
      {
         getfilename(dum,objfile); 
         sprintf(outputfile,"C%s",dum);
      }

      printf("Output filename = %s ...\n", outputfile);

      sub_im=read_nifti_image(objfile, &sub_hdr);
      if( sub_im == NULL)
      {
         printf("Error reading %s, aborting ...\n", objfile);
         exit(0);
      }

      if(!opt_w2d)
      {
         if(sub_hdr.datatype==DT_SIGNED_SHORT || sub_hdr.datatype==DT_UINT16)
         {
            op_im=(char *)computeReslicedImage( (short *)sub_im, sub_hdr, trg_hdr, Xwarp, Ywarp, Zwarp, invT);
         }
         if(sub_hdr.datatype==DT_FLOAT) 
         {
            op_im=(char *)computeReslicedImage( (float *)sub_im, sub_hdr, trg_hdr, Xwarp,Ywarp,Zwarp,invT);
         }
         if(sub_hdr.datatype==DT_UNSIGNED_CHAR) 
         {
            op_im=(char *)computeReslicedImage( (unsigned char *)sub_im, sub_hdr, trg_hdr, Xwarp,Ywarp,Zwarp,invT);
         }
      }

      if(opt_w2d)
      {
         nifti_1_header hdr;
         int nx2d, ny2d, nz2d;
         int np2d, nv2d;
         float dx2d, dy2d, dz2d;
         float *Xw2d, *Yw2d, *Zw2d;
         short *sdum;

         printf("\nReading displacement field: %s\n",w2dfile);

         fp=fopen(w2dfile,"r");
         if(fp==NULL) file_open_error(w2dfile);
         fread(&hdr,sizeof(nifti_1_header),1,fp);
         fread(&extender, sizeof(nifti1_extender), 1, fp);
         if(hdr.dim[0]<1 || hdr.dim[0]>7)
         {
            swapniftiheader(&hdr);
         }
         nx2d = hdr.dim[1];
         ny2d = hdr.dim[2];
         nz2d = hdr.dim[3];
         dx2d = hdr.pixdim[1];
         dy2d = hdr.pixdim[2];
         dz2d = hdr.pixdim[3];

         np2d = nx2d * ny2d;
         nv2d = nx2d * ny2d * nz2d;

         sdum = (short *)calloc(nv2d, sizeof(short));

         Xw2d = (float *)calloc(nv2d ,sizeof(float));
         Yw2d = (float *)calloc(nv2d ,sizeof(float));
         Zw2d = (float *)calloc(nv2d ,sizeof(float));

         fread(sdum,sizeof(short),nv2d,fp);
         for(int i=0; i<nv2d; i++) Xw2d[i] = sdum[i]*hdr.scl_slope;

         fread(sdum,sizeof(short),nv2d,fp);
         for(int i=0; i<nv2d; i++) Yw2d[i] = sdum[i]*hdr.scl_slope;

         fread(sdum,sizeof(short),nv2d,fp);
         for(int i=0; i<nv2d; i++) Zw2d[i] = sdum[i]*hdr.scl_slope;

         delete sdum;
         fclose(fp);

         ///////////////////////////////////////////

         if(sub_hdr.datatype==DT_SIGNED_SHORT || sub_hdr.datatype==DT_UINT16) 
            op_im=(char *)computeReslicedImage((short *)sub_im,sub_hdr,trg_hdr,hdr,Xwarp,Ywarp,Zwarp,invT,Xw2d,Yw2d,Zw2d);

         if(sub_hdr.datatype==DT_FLOAT) 
            op_im=(char *)computeReslicedImage((float *)sub_im,sub_hdr,trg_hdr,hdr,Xwarp,Ywarp,Zwarp,invT,Xw2d,Yw2d,Zw2d);

         if(sub_hdr.datatype==DT_UNSIGNED_CHAR) 
            op_im=(char *)computeReslicedImage((unsigned char *)sub_im,sub_hdr,trg_hdr,hdr,Xwarp,Ywarp,Zwarp,invT,Xw2d,Yw2d,Zw2d);

         free(Xw2d);
         free(Yw2d);
         free(Zw2d);
      }

      trg_hdr.scl_slope = sub_hdr.scl_slope;
      trg_hdr.scl_inter = sub_hdr.scl_inter;
      trg_hdr.datatype = sub_hdr.datatype;

      if(sub_hdr.datatype==DT_SIGNED_SHORT || sub_hdr.datatype==DT_UINT16)
      {
         trg_hdr.bitpix=sizeof(short)*8;
         trg_hdr.dim[0]=3;
         save_nifti_image(outputfile, (short*)op_im, &trg_hdr);
      }
      if(sub_hdr.datatype==DT_FLOAT) 
      {
         trg_hdr.bitpix=sizeof(float)*8;
         trg_hdr.dim[0]=3;
         save_nifti_image(outputfile, (float*)op_im, &trg_hdr);
      }
      if(sub_hdr.datatype==DT_UNSIGNED_CHAR) 
      {
         trg_hdr.bitpix=sizeof(unsigned char)*8;
         trg_hdr.dim[0]=3;
         save_nifti_image(outputfile, (unsigned char*)op_im, &trg_hdr);
      }

      free(op_im);
      free(sub_im);
   }

   free(Xwarp);
   free(Ywarp);
   free(Zwarp);
   free(invT);

   printf("\nEND\n");
}

void print_help_and_exit()
{
	printf("\nUsage: applywarp3d  [ -version -V -h -nn -o <prefix> -T <linear transformation> -w2d <2D displacement field>]\n" 
	"-w <3D displacement field> <input image files ...>\n\n"

	"\tOptional arguments:\n\n"
	"\t-version or -V: Prints software version.\n\n"
	"\t-h or -help: Prints help message.\n\n"
	"\t-nn:  Uses nearest neighbor interpolation (default is linear interpolation).\n\n"
	"\t-o <prefix>: Specifies a prefix for naming output images.\n\n"
	"\t-T <linear transformation>: Specifies a linear transformation matrix be be applied to the input images.\n\n"
	"\t-w2d <2D displacement field>: Specifies a 2D displacement field to be applied to the input images.\n\n"

	"Required arguments:\n\n"
	"\t-w <3D displacement field>: Specifies a 3D displacement field to be applied to the input images.\n\n"
	"\t<input image files ...>: Any number of input image files.\n\n"

	"Notes:\n\n"
	"\tThe input images must in NIFTI format of type 'short' or 'float' or 'unsigned char'.\n\n"
	"\tThe <3D displacement field> is that found by the 3dwarper program.\n\n"
	"\tThe <2D displacement field> can be that found by the unwarp2d program.\n\n"
	"\tIf no <prefix> is specified, the input images will be prepended by the letter 'C' name\n"
	"\tthe output files.\n\n"
);

	exit(0);
}

unsigned char *computeReslicedImage(unsigned char *im1, nifti_1_header hdr1, nifti_1_header hdr2, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;

   float  x,y,z;   
   float  xx,yy,zz;   
   int q;
   int np1;
   unsigned char *im2;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float *beta=NULL, del;
   float *c=NULL;

   nx1 = hdr1.dim[1]; ny1 = hdr1.dim[2]; nz1 = hdr1.dim[3];
   nx2 = hdr2.dim[1]; ny2 = hdr2.dim[2]; nz2 = hdr2.dim[3];
   dx1 = hdr1.pixdim[1]; dy1 = hdr1.pixdim[2]; dz1 = hdr1.pixdim[3];
   dx2 = hdr2.pixdim[1]; dy2 = hdr2.pixdim[2]; dz2 = hdr2.pixdim[3];

	if(opt_cubicspline) 
	{
		beta=computeBeta(&del);
		c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
		cubicSplineAnalysis(im1, c, nx1, ny1, nz1);
	}

	np1=nx1*ny1;

	im2=(unsigned char *)calloc(nx2*ny2*nz2,sizeof(unsigned char));

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

	q=0;
	for(int k=0;k<nz2;k++) 
	{
		for(int j=0;j<ny2;j++) 
		{
  			for(int i=0;i<nx2;i++) 
			{
				zz = k*dz2 - zc2 + Zwarp[q];
				yy = j*dy2 - yc2 + Ywarp[q];
				xx = i*dx2 - xc2 + Xwarp[q];

				x = ( T[0]*xx +T[1]*yy +T[2]*zz  +T[3]   + xc1 )/dx1;
				y = ( T[4]*xx +T[5]*yy +T[6]*zz  +T[7]   + yc1 )/dy1;
				z = ( T[8]*xx +T[9]*yy +T[10]*zz +T[11]  + zc1 )/dz1;

				if(opt_cubicspline)
					im2[q++] = (unsigned char) cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del);
				else if(opt_nn)
	   				im2[q++] = (uchar) ( nearestNeighbor(x,y,z,im1,nx1,ny1,nz1,np1) + 0.5);
				else
	   				im2[q++] = (uchar) ( linearInterpolator(x,y,z,im1,nx1,ny1,nz1,np1) + 0.5);
			}
		}
	}

	if(opt_cubicspline) 
	{
		free(beta);
		free(c);
	}

	return( im2 );
}

short *computeReslicedImage(short *im1, nifti_1_header hdr1, nifti_1_header hdr2, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;

   float  x,y,z;   
   float  xx,yy,zz;   
   int q;
   int np1;
   short *im2;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float *beta=NULL, del;
   float *c=NULL;

   nx1 = hdr1.dim[1]; ny1 = hdr1.dim[2]; nz1 = hdr1.dim[3];
   nx2 = hdr2.dim[1]; ny2 = hdr2.dim[2]; nz2 = hdr2.dim[3];
   dx1 = hdr1.pixdim[1]; dy1 = hdr1.pixdim[2]; dz1 = hdr1.pixdim[3];
   dx2 = hdr2.pixdim[1]; dy2 = hdr2.pixdim[2]; dz2 = hdr2.pixdim[3];

   //printf("Target voxel size: %f %f %f\n",dx2,dy2,dz2);
   //printf("Target matrix size: %d %d %d\n",nx2,ny2,nz2);

   //printf("Subject voxel size: %f %f %f\n",dx1,dy1,dz1);
   //printf("Subject matrix size: %d %d %d\n",nx1,ny1,nz1);

   if(opt_cubicspline) 
   {
      beta=computeBeta(&del);
      c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
      cubicSplineAnalysis(im1, c, nx1, ny1, nz1);
   }

   np1=nx1*ny1;

   im2=(short *)calloc(nx2*ny2*nz2,sizeof(short));
   if( im2 == NULL )
   {
      memory_allocation_error("im2");
   }

   xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
   yc1=dy1*(ny1-1)/2.0;
   zc1=dz1*(nz1-1)/2.0;

   xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
   yc2=dy2*(ny2-1)/2.0;
   zc2=dz2*(nz2-1)/2.0;

   q=0;
   for(int k=0;k<nz2;k++) 
   {
      for(int j=0;j<ny2;j++) 
      {
         for(int i=0;i<nx2;i++) 
         {
            xx = i*dx2 - xc2 + Xwarp[q];
            yy = j*dy2 - yc2 + Ywarp[q];
            zz = k*dz2 - zc2 + Zwarp[q];

            x = ( T[0]*xx +T[1]*yy +T[2]*zz  +T[3]   + xc1 )/dx1;
            y = ( T[4]*xx +T[5]*yy +T[6]*zz  +T[7]   + yc1 )/dy1;
            z = ( T[8]*xx +T[9]*yy +T[10]*zz +T[11]  + zc1 )/dz1;

            if(opt_cubicspline)
               im2[q++] = (short) (cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
            else if(opt_nn)
               im2[q++] = (short) (nearestNeighbor(x,y,z,im1,nx1,ny1,nz1,np1)+0.5);
            else
               im2[q++] = (short) (linearInterpolator(x,y,z,im1,nx1,ny1,nz1,np1)+0.5);
         }
      }
   }

//printf("%d %d %d\n", q, nx2*ny2*nz2, im2[(nx2*ny2)*nz2/2 + nx2*ny2/2 + nx2/2]);

   if(opt_cubicspline) 
   {
      free(beta);
      free(c);
   }

   return( im2 );
}

float *computeReslicedImage(float *im1, nifti_1_header hdr1, nifti_1_header hdr2, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;

   float  x,y,z;   
   float  xx,yy,zz;   
   int q;
   int np1;
   float *im2;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float *beta=NULL, del;
   float *c=NULL;

   nx1 = hdr1.dim[1]; ny1 = hdr1.dim[2]; nz1 = hdr1.dim[3];
   nx2 = hdr2.dim[1]; ny2 = hdr2.dim[2]; nz2 = hdr2.dim[3];
   dx1 = hdr1.pixdim[1]; dy1 = hdr1.pixdim[2]; dz1 = hdr1.pixdim[3];
   dx2 = hdr2.pixdim[1]; dy2 = hdr2.pixdim[2]; dz2 = hdr2.pixdim[3];

	if(opt_cubicspline) 
	{
		beta=computeBeta(&del);
		c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
		cubicSplineAnalysis(im1, c, nx1, ny1, nz1);
	}

   np1=nx1*ny1;

   im2=(float *)calloc(nx2*ny2*nz2,sizeof(float));

   xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
   yc1=dy1*(ny1-1)/2.0;
   zc1=dz1*(nz1-1)/2.0;

   xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
   yc2=dy2*(ny2-1)/2.0;
   zc2=dz2*(nz2-1)/2.0;

	q=0;
	for(int k=0;k<nz2;k++) 
	{
		for(int j=0;j<ny2;j++) 
		{
  			for(int i=0;i<nx2;i++) 
			{
				zz = k*dz2 - zc2 + Zwarp[q];
				yy = j*dy2 - yc2 + Ywarp[q];
				xx = i*dx2 - xc2 + Xwarp[q];

				x = ( T[0]*xx +T[1]*yy +T[2]*zz  +T[3]   + xc1 )/dx1;
				y = ( T[4]*xx +T[5]*yy +T[6]*zz  +T[7]   + yc1 )/dy1;
				z = ( T[8]*xx +T[9]*yy +T[10]*zz +T[11]  + zc1 )/dz1;

				if(opt_cubicspline)
					im2[q++] = cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del);
				else if(opt_nn)
	   				im2[q++]=nearestNeighbor(x,y,z,im1,nx1,ny1,nz1,np1);
				else
	   				im2[q++]=linearInterpolator(x,y,z,im1,nx1,ny1,nz1,np1);
			}
		}
	}

	if(opt_cubicspline) 
	{
		free(beta);
		free(c);
	}

   return( im2 );
}

short *computeReslicedImage(short *im1, nifti_1_header hdr1, nifti_1_header hdr2, nifti_1_header hdr3,
float *Xwarp, float *Ywarp, float *Zwarp, float *T, float *Xw2d, float *Yw2d, float *Zw2d)
{
   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   int nx3, ny3, nz3;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;
   float dx3, dy3, dz3;

   float  x,y,z;   
   float  xx,yy,zz;   
   int q;
   int np1,np2,np3;
   short *im2;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float xc3, yc3, zc3;
   float wx,wy,wz;

	float *beta=NULL, del;
	float *c1=NULL,*c2=NULL,*c3=NULL, *c4=NULL;

   nx1 = hdr1.dim[1]; ny1 = hdr1.dim[2]; nz1 = hdr1.dim[3];
   nx2 = hdr2.dim[1]; ny2 = hdr2.dim[2]; nz2 = hdr2.dim[3];
   nx3 = hdr3.dim[1]; ny3 = hdr3.dim[2]; nz3 = hdr3.dim[3];
   dx1 = hdr1.pixdim[1]; dy1 = hdr1.pixdim[2]; dz1 = hdr1.pixdim[3];
   dx2 = hdr2.pixdim[1]; dy2 = hdr2.pixdim[2]; dz2 = hdr2.pixdim[3];
   dx3 = hdr3.pixdim[1]; dy3 = hdr3.pixdim[2]; dz3 = hdr3.pixdim[3];

	if(opt_cubicspline) 
	{
		beta=computeBeta(&del);
		c1 = (float *)calloc(nx3*ny3*nz3, sizeof(float));
		c2 = (float *)calloc(nx3*ny3*nz3, sizeof(float));
		c3 = (float *)calloc(nx3*ny3*nz3, sizeof(float));
		c4 = (float *)calloc(nx1*ny1*nz1, sizeof(float));
		cubicSplineAnalysis(Xw2d, c1, nx3, ny3, nz3);
		cubicSplineAnalysis(Yw2d, c2, nx3, ny3, nz3);
		cubicSplineAnalysis(Zw2d, c3, nx3, ny3, nz3);
		cubicSplineAnalysis(im1,  c4, nx1, ny1, nz1);
	}

	np1=nx1*ny1;
	np2=nx2*ny2;
	np3=nx3*ny3;

	im2=(short *)calloc(nx2*ny2*nz2,sizeof(short));

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

   	xc3=dx3*(nx3-1)/2.0;     /* +---+---+ */
	yc3=dy3*(ny3-1)/2.0;
	zc3=dz3*(nz3-1)/2.0;

	q=0;
	for(int k=0;k<nz2;k++) 
	{
		for(int j=0;j<ny2;j++) 
		{
  			for(int i=0;i<nx2;i++) 
			{
				xx = i*dx2 - xc2 + Xwarp[q];
				yy = j*dy2 - yc2 + Ywarp[q];
				zz = k*dz2 - zc2 + Zwarp[q];

				x = ( T[0]*xx +T[1]*yy +T[2]*zz  +T[3] );
				y = ( T[4]*xx +T[5]*yy +T[6]*zz  +T[7] );
				z = ( T[8]*xx +T[9]*yy +T[10]*zz +T[11]);

				xx = (x  + xc3)/dx3;
				yy = (y  + yc3)/dy3;
				zz = (z  + zc3)/dz3;

				if(opt_cubicspline)
				{
					wx = cubicSplineSynthesis(c1, nx3, ny3, nz3, xx, yy, zz, beta, del);
					wy = cubicSplineSynthesis(c2, nx3, ny3, nz3, xx, yy, zz, beta, del);
					wz = cubicSplineSynthesis(c3, nx3, ny3, nz3, xx, yy, zz, beta, del);
				}
				else
				{
					wx = linearInterpolator(xx,yy,zz,Xw2d,nx3,ny3,nz3,np3);
					wy = linearInterpolator(xx,yy,zz,Yw2d,nx3,ny3,nz3,np3);
					wz = linearInterpolator(xx,yy,zz,Zw2d,nx3,ny3,nz3,np3);
				}

				xx = ( x + wx + xc1 )/dx1;
				yy = ( y + wy + yc1 )/dy1;
				zz = ( z + wz + zc1 )/dz1;

				if(opt_cubicspline)
					im2[q++] = (short)(cubicSplineSynthesis(c4, nx1, ny1, nz1, xx, yy, zz, beta, del)+0.5);
				else if(opt_nn)
	   				im2[q++] = (short)(nearestNeighbor(xx,yy,zz,im1,nx1,ny1,nz1,np1)+0.5);
				else
	   				im2[q++] = (short)(linearInterpolator(xx,yy,zz,im1,nx1,ny1,nz1,np1)+0.5);
			}
		}
	}

	if(opt_cubicspline) 
	{
		free(beta);
		free(c1);
		free(c2);
		free(c3);
		free(c4);
	}

	return( im2 );
}

float *computeReslicedImage(float *im1, nifti_1_header hdr1, nifti_1_header hdr2, nifti_1_header hdr3,
float *Xwarp, float *Ywarp, float *Zwarp, float *T, float *Xw2d, float *Yw2d, float *Zw2d)
{
   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   int nx3, ny3, nz3;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;
   float dx3, dy3, dz3;

  	float  x,y,z;   
  	float  xx,yy,zz;   
	int q;
	int np1,np2,np3;
	float *im2;
	float xc1, yc1, zc1;
	float xc2, yc2, zc2;
	float xc3, yc3, zc3;
	float wx,wy,wz;

   nx1 = hdr1.dim[1]; ny1 = hdr1.dim[2]; nz1 = hdr1.dim[3];
   nx2 = hdr2.dim[1]; ny2 = hdr2.dim[2]; nz2 = hdr2.dim[3];
   nx3 = hdr3.dim[1]; ny3 = hdr3.dim[2]; nz3 = hdr3.dim[3];
   dx1 = hdr1.pixdim[1]; dy1 = hdr1.pixdim[2]; dz1 = hdr1.pixdim[3];
   dx2 = hdr2.pixdim[1]; dy2 = hdr2.pixdim[2]; dz2 = hdr2.pixdim[3];
   dx3 = hdr3.pixdim[1]; dy3 = hdr3.pixdim[2]; dz3 = hdr3.pixdim[3];

	np1=nx1*ny1;
	np2=nx2*ny2;
	np3=nx3*ny3;

	im2=(float *)calloc(nx2*ny2*nz2,sizeof(float));

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

   	xc3=dx3*(nx3-1)/2.0;     /* +---+---+ */
	yc3=dy3*(ny3-1)/2.0;
	zc3=dz3*(nz3-1)/2.0;

	q=0;
	for(int k=0;k<nz2;k++) 
	{
		for(int j=0;j<ny2;j++) 
		{
  			for(int i=0;i<nx2;i++) 
			{
				xx = i*dx2 - xc2 + Xwarp[q];
				yy = j*dy2 - yc2 + Ywarp[q];
				zz = k*dz2 - zc2 + Zwarp[q];

				x = ( T[0]*xx +T[1]*yy +T[2]*zz  +T[3] );
				y = ( T[4]*xx +T[5]*yy +T[6]*zz  +T[7] );
				z = ( T[8]*xx +T[9]*yy +T[10]*zz +T[11]);

				xx = (x  + xc3)/dx3 ;
				yy = (y  + yc3)/dy3 ;
				zz = (z  + zc3)/dz3 ;

				wx = linearInterpolator(xx,yy,zz,Xw2d,nx3,ny3,nz3,np3);
				wy = linearInterpolator(xx,yy,zz,Yw2d,nx3,ny3,nz3,np3);
				wz = linearInterpolator(xx,yy,zz,Zw2d,nx3,ny3,nz3,np3);

				xx = ( x + wx + xc1 )/dx1;
				yy = ( y + wy + yc1 )/dy1;
				zz = ( z + wz + zc1 )/dz1;

				if(opt_nn)
	   				im2[q++]=nearestNeighbor(xx,yy,zz,im1,nx1,ny1,nz1,np1);
				else
	   				im2[q++]=linearInterpolator(xx,yy,zz,im1,nx1,ny1,nz1,np1);
			}
		}
	}

	return( im2 );
}

unsigned char *computeReslicedImage(unsigned char *im1, nifti_1_header hdr1, nifti_1_header hdr2, nifti_1_header hdr3,
float *Xwarp, float *Ywarp, float *Zwarp, float *T, float *Xw2d, float *Yw2d, float *Zw2d)
{
   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   int nx3, ny3, nz3;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;
   float dx3, dy3, dz3;

   float  x,y,z;   
   float  xx,yy,zz;   
   int q;
   int np1,np2,np3;
   unsigned char *im2;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float xc3, yc3, zc3;
   float wx,wy,wz;

   nx1 = hdr1.dim[1]; ny1 = hdr1.dim[2]; nz1 = hdr1.dim[3];
   nx2 = hdr2.dim[1]; ny2 = hdr2.dim[2]; nz2 = hdr2.dim[3];
   nx3 = hdr3.dim[1]; ny3 = hdr3.dim[2]; nz3 = hdr3.dim[3];
   dx1 = hdr1.pixdim[1]; dy1 = hdr1.pixdim[2]; dz1 = hdr1.pixdim[3];
   dx2 = hdr2.pixdim[1]; dy2 = hdr2.pixdim[2]; dz2 = hdr2.pixdim[3];
   dx3 = hdr3.pixdim[1]; dy3 = hdr3.pixdim[2]; dz3 = hdr3.pixdim[3];

	np1=nx1*ny1;
	np2=nx2*ny2;
	np3=nx3*ny3;

	im2=(uchar *)calloc(nx2*ny2*nz2,sizeof(unsigned char));

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

   	xc3=dx3*(nx3-1)/2.0;     /* +---+---+ */
	yc3=dy3*(ny3-1)/2.0;
	zc3=dz3*(nz3-1)/2.0;

	q=0;
	for(int k=0;k<nz2;k++) 
	{
		for(int j=0;j<ny2;j++) 
		{
  			for(int i=0;i<nx2;i++) 
			{
				xx = i*dx2 - xc2 + Xwarp[q];
				yy = j*dy2 - yc2 + Ywarp[q];
				zz = k*dz2 - zc2 + Zwarp[q];

				x = ( T[0]*xx +T[1]*yy +T[2]*zz  +T[3] );
				y = ( T[4]*xx +T[5]*yy +T[6]*zz  +T[7] );
				z = ( T[8]*xx +T[9]*yy +T[10]*zz +T[11]);

				xx = (x  + xc3)/dx3 ;
				yy = (y  + yc3)/dy3 ;
				zz = (z  + zc3)/dz3 ;

				wx = linearInterpolator(xx,yy,zz,Xw2d,nx3,ny3,nz3,np3);
				wy = linearInterpolator(xx,yy,zz,Yw2d,nx3,ny3,nz3,np3);
				wz = linearInterpolator(xx,yy,zz,Zw2d,nx3,ny3,nz3,np3);

				xx = ( x + wx + xc1 )/dx1;
				yy = ( y + wy + yc1 )/dy1;
				zz = ( z + wz + zc1 )/dz1;

				if(opt_nn)
	   				im2[q++]=(uchar)(nearestNeighbor(xx,yy,zz,im1,nx1,ny1,nz1,np1)+0.5);
				else
	   				im2[q++]=(uchar)(linearInterpolator(xx,yy,zz,im1,nx1,ny1,nz1,np1)+0.5);
			}
		}
	}

	return( im2 );
}
