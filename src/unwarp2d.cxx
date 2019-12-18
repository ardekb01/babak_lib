// Revision 5/13/11: Smoothing functions (smoothX, smoothY, etc.) were transferred to smooth.cxx under babak_lib

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
#include <smooth.h>
#include <minmax.h>

#define YES 1
#define NO 0

extern float *resizeXY(float *image1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2);
extern short *resizeXY(short *image1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2);

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-v", 0,  'v'},
   {"-sub", 1,  's'},
   {"-obj", 1,  's'},

   {"-I", 0,  'I'},
   {"-T", 1,  'T'},

   {"-o", 1,  'o'},

   {"-sd", 1,  'S'},
   {"-trg", 1,  't'},
   {"-iter", 1,  '5'},
   {"-u", 1,  '6'},
   {"-w", 1,  '7'},
   {"-s", 1,  '9'},
   {0, 0,  0}
};

int opt_I=NO;
int opt_sd=NO;
int opt_sub=NO;
int opt_trg=NO;
int opt_w=NO;
int opt_iter=NO;
int opt_s=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit();
void check_W_permission(char *file);
void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz);
void swapN(char *in, int N);
void swapByteOrder(char *in, int N);

int main(int argc, char **argv)
{
  nifti_1_header trg_hdr;
  nifti_1_header sub_hdr;

  float T[16];

	FILE *fp;

	float Vt,Vo;
	float Sx, Sx2, Sy, Sy2, Sxy;
	float num,den;

	float sd=1;
	int xopt, yopt;
	int Wx,Wy;
	int Lx,Ly;
	int N=0;	// N=2*L+1
	int N2=0;	// N*N
	int niter=0;

	int Onx,Ony,Onz,Onp;
	int Tnx,Tny,Tnz,Tnp,Tnv;
	int HRnx, HRny;

	float Odx,Ody,Odz;
	float Tdx,Tdy,Tdz;
	float HRdx, HRdy;

	float *Xw, *Yw;
	float *Xwarp, *Ywarp, *Zwarp, *Xwarp_ptr, *Ywarp_ptr;
	float *Xtmp, *Ytmp;
	float CC, CCMAX;

	char *dum;
   char outputfile[1024]="";
	char trgfile[1024];

   char subfile[1024];
   char warpfile[1024];

   char Tfile[1024]="";

	short *trg;
	short *objOrig, *obj;

	short *HRtrg;	// high res. target and object images
	short *HRobj;

	float *ARobj;	// array extracted from object and target images
	float *ARtrg;

	int window_width=1;

   dum = (char *)malloc(1024);
   warpfile[0]='\0';
   subfile[0]='\0';

	while( (opt=getoption(argc, argv, options)) != -1)
	{
      switch (opt) {
         case 'I':
            opt_I=YES;
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'T':
            sprintf(Tfile,"%s",optarg);
            break;
         case 's':
            sprintf(subfile,"%s",optarg);
            opt_sub=YES;
            break;
         case 't':
            sprintf(trgfile,"%s",optarg);
            opt_trg=YES;
            break;
         case 'o':
            sprintf(outputfile,"%s",optarg);
            break;
			case 'S':
				sd=atof(optarg);
				opt_sd=YES;
				break;
			case '5':
				niter=atoi(optarg);
				opt_iter=YES;
				break;
			case '6':
				sprintf(warpfile,"%s",optarg);
				break;
			case '7':
				N=atoi(optarg);
				opt_w=YES;
				break;
			case '9':
				window_width=atoi(optarg);
				opt_s=YES;
				break;
			case '?':
				print_help_and_exit();
		}
	}


	if( !opt_sub || !opt_trg )
		print_help_and_exit();

	if(!opt_w || N<=0) N=5;
	if( (N%2)==0 ) N+=1;			// make sure it's odd
	
	if(!opt_s || window_width<=0) window_width=5;
	if( (window_width%2)==0 ) window_width+=1;			// make sure it's odd

	if(!opt_iter || niter<=0 ) niter=4;
	printf("\nNumber of iterations = %d\n",niter);
	
   //////////////////////////////////////////////////////////////////////////////////////////
   /// If no output filename is specified, the program automatically sets it by putting the
   /// letter 'C' in front of the object filename. For example, if the object filename is 
   /// test.nii, the output filename will be Ctest.nii.
   //////////////////////////////////////////////////////////////////////////////////////////
   if(outputfile[0]=='\0') 
   {
      // returns the filename in dum, without the path information.
      getfilename(dum,subfile); 
      sprintf(outputfile,"C%s",dum);
   }

   if(opt_v)
   {
      printf("Output (unwarped) image file = %s\n",outputfile);
   }

   /////////////////////////////////////////////////////////////////////////////////////////

   if(opt_v)
   {
      printf("Target image file = %s\n",trgfile);
   }
   trg=(short *)read_nifti_image(trgfile, &trg_hdr);
   Tnx=trg_hdr.dim[1]; Tny=trg_hdr.dim[2]; Tnz=trg_hdr.dim[3];
   Tdx=trg_hdr.pixdim[1]; Tdy=trg_hdr.pixdim[2]; Tdz=trg_hdr.pixdim[3];

   Tnv=Tnx*Tny*Tnz;
   Tnp=Tnx*Tny;

   if(trg_hdr.datatype != DT_SIGNED_SHORT && trg_hdr.datatype != 512)
   {
      printf("\nSorry, this program currently only handles images of type DT_SIGNED_SHORT\n.");
      exit(0);
   }
   /////////////////////////////////////////////////////////////////////////////////////////

   if(opt_v)
   {
      printf("Subject image file = %s\n",subfile);
   }
   objOrig=(short *)read_nifti_image(subfile, &sub_hdr);
   Onx=sub_hdr.dim[1]; Ony=sub_hdr.dim[2]; Onz=sub_hdr.dim[3];
   Odx=sub_hdr.pixdim[1]; Ody=sub_hdr.pixdim[2]; Odz=sub_hdr.pixdim[3];
   Onp=Onx*Ony;

   if(sub_hdr.datatype != DT_SIGNED_SHORT && sub_hdr.datatype != 512)
   {
      printf("\nSorry, this program currently only handles images of type DT_SIGNED_SHORT\n.");
      exit(0);
   }

   /////////////////////////////////////////////////////////////////////////////////////////

   for(int i=0; i<16; i++) T[i]=0.0;
   T[0]=T[5]=T[10]=T[15]=1.0;
   if( Tfile[0]!='\0' ) 
   {
      if(opt_v)
      {
         printf("Reading subject->target rigid-body transformation from file: %s ...\n", Tfile);
      }

      loadTransformation(Tfile, T);
   }
   else if(!opt_I)
   {
      if(opt_v)
      {
         printf("Computing subject->target rigid-body transformation from header information ...\n");
      }

      sub2trg_rigid_body_transformation(T, subfile, trgfile);
   }

   if(opt_v)
   {
      printMatrix(T,4,4,"Subject->target rigid-body transformation:",NULL);
   }

   /////////////////////////////////////////////////////////////////////////////////////////
   {
      float *invT;

      invT=inv4(T);
      obj=resliceImage(objOrig, Onx, Ony, Onz, Odx, Ody, Odz, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, invT, LIN);
      free(invT);
   }

   /////////////////////////////////////////////////////////////////////////////////////////

	Lx = (N-1)/2;
	Ly = (int)(Lx*Tdx/Tdy + 0.5);
	printf("\nCorrelation window size = %d x %d\n",2*Lx+1, 2*Ly+1);
	N2 = (2*Lx+1)*(2*Ly+1);

	Wx = (window_width-1)/2;
	Wy = (int)(Wx*Tdx/Tdy + 0.5);
	printf("\nInitial search window = %d x %d \n",2*Wx+1, 2*Wy+1);

	ARobj=(float *)calloc(N2,sizeof(float));
	ARtrg=(float *)calloc(N2,sizeof(float));

	if(!sd || sd<0.0) sd = 2.0*Wx;

	Xwarp=(float *)calloc(Tnv,sizeof(float));
	Ywarp=(float *)calloc(Tnv,sizeof(float));
	Zwarp=(float *)calloc(Tnv,sizeof(float));

	for(int s=0; s<Tnz; s++)
	{
		Xwarp_ptr = Xwarp+s*Tnp;
		Ywarp_ptr = Ywarp+s*Tnp;

		printf("\nProcessing slice %d\n",s+1);

		for(int n=0; n<Tnp; n++) Xwarp_ptr[n]=Ywarp_ptr[n]=0.0;

		for(int iter=0; iter<niter; iter++)
		{
			printf("\n\tIteration %d\n",iter);

			HRdx = (float)(Tdx*pow(2.0,(niter-iter-1.0)));
			HRdy = (float)(Tdy*pow(2.0,(niter-iter-1.0)));

			HRnx = (int)(Tnx/pow(2.0,(niter-iter-1.0)) + 0.5 );
			HRny = (int)(Tny/pow(2.0,(niter-iter-1.0)) + 0.5 );


			printf("\tMatrix size = %d x %d (voxels)\n", HRnx, HRny);
			printf("\tVoxel size = %8.6f x %8.6f (mm3)\n", HRdx,HRdy);

			Xw=resizeXY(Xwarp_ptr, Tnx, Tny, Tdx, Tdy, HRnx, HRny, HRdx, HRdy);
			Yw=resizeXY(Ywarp_ptr, Tnx, Tny, Tdx, Tdy, HRnx, HRny, HRdx, HRdy);

			{
				float *tmp;
				float StandDev;

				StandDev = (0.5/log(2.0)) * ( HRdx*HRdx - Odx*Odx );

				if(StandDev>0.0)
				{
                    if( ( HRdx*HRdx - Tdx*Tdx ) > 0.0 )
                    {
					   StandDev=(float)( sqrt( (0.5/log(2.0)) * ( HRdx*HRdx - Tdx*Tdx ) )/Tdx );
                    }

					tmp = smoothXY(obj+s*Tnp,Tnx,Tny,StandDev);
					HRobj=computeReslicedImage(tmp, Tnx, Tny, Tdx, Tdy, HRnx, HRny, HRdx, HRdy, Xw, Yw);
					free(tmp);
				}
				else
				{
					HRobj=computeReslicedImage(obj+s*Tnp, Tnx, Tny, Tdx, Tdy, HRnx, HRny, HRdx, HRdy, Xw, Yw);
				}
			}

			HRtrg=resizeXY(trg+s*Tnp, Tnx ,Tny, Tdx, Tdy, HRnx, HRny, HRdx, HRdy);

			for(int j=0; j<HRny; j++)
			for(int i=0; i<HRnx; i++)
			{
				extractArray(HRtrg, HRnx, HRny, i, j, Lx, Ly, ARtrg);

				Sy=Sy2=0.0;
				for(int n=0; n<N2; n++)
				{
					Vt = ARtrg[n];
					Sy += Vt;
					Sy2 += (Vt*Vt);
				}
	
				if( Sy == 0.0 )
				{
					Xw[j*HRnx + i] = 0.0;
					Yw[j*HRnx + i] = 0.0;
					continue;
				}
	
				CCMAX=0.0; 	// IMPORTANT: we are not interested in -tive correlations
							// if CMAX is set to -1, program given unexpected results
	
				xopt=yopt=0;
	
				for(int x=-Wx; x<=Wx; x++)
				for(int y=-Wy; y<=Wy; y++)
				{
					extractArray(HRobj, HRnx, HRny, i+x, j+y, Lx, Ly, ARobj);
	
					Sx=Sx2=Sxy=0.0;
					for(int n=0; n<N2; n++)
					{
						Vo = ARobj[n];
						Sx += Vo;
						Sx2 += (Vo*Vo);
						Sxy += (Vo*ARtrg[n]);
					}
	
					num = Sxy-Sx*Sy/N2;
					den = (float)sqrt( (double)(Sx2-Sx*Sx/N2)*(Sy2-Sy*Sy/N2) );
	
					if(den==0.0) continue;
	
					CC = num/den;
	
					if( CC>CCMAX ) { CCMAX=CC; xopt=x; yopt=y; }
				}
	
				Xw[j*HRnx + i] = xopt*HRdx;
				Yw[j*HRnx + i] = yopt*HRdy;
			}
	
			Xtmp=smoothXY(Xw, HRnx, HRny, sd);
			Ytmp=smoothXY(Yw, HRnx, HRny, sd);

			free(Xw); free(Yw);

			Xw=resizeXY(Xtmp, HRnx, HRny, HRdx, HRdy, Tnx, Tny, Tdx, Tdy);
			Yw=resizeXY(Ytmp, HRnx, HRny, HRdx, HRdy, Tnx, Tny, Tdx, Tdy);
	
			free(Xtmp); free(Ytmp);
	
			for(int n=0; n<Tnp; n++)
			{
				Xwarp_ptr[n] += Xw[n];
				Ywarp_ptr[n] += Yw[n];
			}
	
			free(Xw); free(Yw);

			free(HRobj);
			free(HRtrg);
		}
	}

   ///////////////////////////////////////////////////////////////////////////////////////////////////////
   {
      short *o;
      float min=0.0, max=0.0, s=0.0;
      short *sdum;
      nifti1_extender extender;

      //////////////////////////////////////////////////////////////////////////////////////////
      /// If no displacement field filename is specified, the program automatically sets it to be
      /// the subject filename appended by '_wrp'. For example, if the subject filename is
      /// test.nii, the warp parameters filename will be test_wrp.nii
      //////////////////////////////////////////////////////////////////////////////////////////
      if(warpfile[0]=='\0')
      {
         // returns the filename in dum, without the path information.
         getfilename(dum,subfile);

         sprintf(warpfile,"%s_wrp.nii", strsep(&dum, "."));
      }

      if(opt_v)
      {
         printf("\nOutput warp parameters file = %s\n",warpfile);
      }

      combine_warps_and_trans(Tnx, Tny, Tnz, Tdx, Tdy, Tdz, Xwarp, Ywarp, Zwarp, T);

      o=computeReslicedImage(objOrig, Onx, Ony, Onz, Odx, Ody, Odz, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, Xwarp, Ywarp, Zwarp);
      save_nifti_image(outputfile, o, &trg_hdr);
      free(o);

      minmax(Xwarp, Tnv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Ywarp, Tnv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Zwarp, Tnv, min, max);
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

      sdum = (short *)calloc(Tnv, sizeof(short));
      for(int i=0; i<Tnv; i++)
      {
         sdum[i] = (short)(Xwarp[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),Tnv,fp);

      for(int i=0; i<Tnv; i++)
      {
         sdum[i] = (short)(Ywarp[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),Tnv,fp);

      for(int i=0; i<Tnv; i++)
      {
         sdum[i] = (short)(Zwarp[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),Tnv,fp);

      fclose(fp);

      free(sdum);
   }
   ///////////////////////////////////////////////////////////////////////////////////////////////////////

   free(Xwarp);
   free(Ywarp);
   free(Zwarp);
}

void print_help_and_exit()
{
	printf("\n\nUsage: unwarp2d -obj <object image> -trg <target image>\n"
	"[-o <output image filename>] [-u <unwarping parameters file>] [-v] [-T <transformation matrix>]\n"
	"[-iter <integer>] [-w <integer>] [-s <integer>]\n\n"
	"-T is used to specify a 4x4 linear transformation matrix. This transformation is applied to the object image\n"
	"first and the resulting image is matched slice-by-slice to the target image.\n\n"
	"-v specifies vertical phase encoding (default=horizontal)\n\n"
	"-o where to save the unwarped image (default=C<object_filename>)\n\n"
	"-u where to save the unwarping parameters (default=object_filename.wrp)\n\n"
	"-w correlation window size (default=31)\n\n"
	"-s search window size (default=21)\n\n"
	"-iter specifies number of iterations (default=5)\n\n");
	exit(0);
}

void check_W_permission(char *file)
{
	if( access(file,F_OK) == 0 && access(file,W_OK) == -1 )
	{
		printf("\n\nWrite permission for %s denied! Aborting ...\n\n",file);
		exit(0);
	}
}

void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
	if ( hdr.hk.sizeof_hdr != 348 )
		swapByteOrder( (char *) &hdr.dime.datatype, sizeof(short) );

	if( hdr.dime.datatype != 4 && hdr.dime.datatype != 512)
	{
		printf("\n\nSorry, but this program can only handle data type 'short'! Aborting ...\n\n");
		exit(0);

	}

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapN( (char *) hdr.dime.dim , 16);

		for(int i=0; i<8; i++)
			swapByteOrder( (char *) &hdr.dime.pixdim[i], sizeof(float) );
	}

	*ny=hdr.dime.dim[2];
	*nx=hdr.dime.dim[1];
	*nz=hdr.dime.dim[3];

	printf("Matrix size = %d x %d x %d (voxels)\n", *nx, *ny, *nz);

	*dx=hdr.dime.pixdim[1];
	*dy=hdr.dime.pixdim[2];
	*dz=hdr.dime.pixdim[3];

        printf("Voxel size = %8.6f x %8.6f x %8.6f (mm3)\n", *dx,*dy,*dz);
}
