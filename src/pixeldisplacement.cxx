
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <unistd.h>
#include <spm_analyze.h>
#include <babak_lib.h>
#include <sph.h>
#include <landmarks.h>
#include <minmax.h>
#include <ctype.h>
#include <stats.h>

#define YES 1
#define NO 0

/////////////////////////////////////
// Global variables required by ATRA
int PILcloudthreshold=0;
/////////////////////////////////////

int opt;

static struct option options[] =
{
   {"-lm", 1, 'l'},
   {"-T", 1, 'T'},
   {"-i", 1, 'i'},
   {"-o",1,'o'},
   {"-v", 0, 'v'},
   {"-h", 0, 'h'},
   {"-help", 0, 'h'},
   {0, 0, 0}
};

float *pixeldisp(DIM dim, float *T)
{
  float xc,yc,zc; 
  float Ax,Bx;
  float Ay,By;
  float Az,Bz;
  float xx,yy,zz; /* translation parameters */
  float x,y,z;   
  float *im;
  int q;

  im=(float *)calloc(dim.nx*dim.ny*dim.nz,sizeof(float));

  xc=dim.dx*(dim.nx-1)/2.0;      /* +---+---+ */
  yc=dim.dy*(dim.ny-1)/2.0;
  zc=dim.dz*(dim.nz-1)/2.0;

  q=0;
  for(int k=0;k<dim.nz;k++) 
  {
    zz=k*dim.dz-zc;
    Bx=T[2]*zz+T[3];
    By=T[6]*zz+T[7];
    Bz=T[10]*zz+T[11];
    for(int j=0;j<dim.ny;j++) 
    {
      yy=j*dim.dy-yc;
      Ax=T[1]*yy+Bx;
      Ay=T[5]*yy+By;
      Az=T[9]*yy+Bz;

      for(int i=0;i<dim.nx;i++) 
      {
        xx=i*dim.dx-xc;

        x=T[0]*xx+Ax;
        y=T[4]*xx+Ay;
        z=T[8]*xx+Az;

        im[q++] = sqrtf( (x-xx)*(x-xx) + (y-yy)*(y-yy) + (z-zz)*(z-zz) );
      }
    }
  }

  return(im);
}

int main(int argc, char **argv)
{
  float T[16];
  short *PILbraincloud=NULL;
  nifti_1_header PILbraincloud_hdr;
  nifti_1_header ipim_hdr;
  DIM PILbraincloud_dim;
  DIM ipim_dim;
  char temporaryFilename[DEFAULT_STRING_LENGTH]; // place holder for temporary filenames
  char ipimfile[DEFAULT_STRING_LENGTH]=""; 
  char transFile[DEFAULT_STRING_LENGTH]=""; 
  char lmfile[DEFAULT_STRING_LENGTH]="";

  short *ipim=NULL;

  opt_png=NO;
  opt_txt=NO;

  while( (opt=getoption(argc, argv, options)) != -1)
  {
    switch (opt) 
    {
      case 'T':
        sprintf(transFile,"%s",optarg);
        break;
      case 'i':
        strcpy(ipimfile,optarg);
        break;
      case 'v':
        opt_v=YES;
        break;
      case 'l':
        strcpy(lmfile,optarg);
        break;
      case '?':
        break;
    }
  }

  getARTHOME();

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //Receive transformation matrix
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if(transFile[0]=='\0')
  {
    printf("Please specify the transformation matrix using the \"-T <matrix>\" flag.\n");
    exit(0);
  }

  if(opt_v)
  {
    printf("Transformation matrix = %s\n",transFile);
  }

  // ensure the input transformatin file exists, is readable, and has the expected size
  if ( checkFileExistence(transFile)==0 )
  {
    printf("Error: File %s does not exist! Aborting ...\n",transFile);
    exit(0);
  }

  if ( checkFileReadOK(transFile)==0 )
  {
    printf("Error: Read permission for %s denied! Aborting ...\n",transFile);
    exit(0);
  }

  loadTransformation(transFile, T);

  if(opt_v)
  {
    printf("%9.6f  %9.6f  %9.6f  %9.6f\n",T[0],T[1],T[2],T[3]);
    printf("%9.6f  %9.6f  %9.6f  %9.6f\n",T[4],T[5],T[6],T[7]);
    printf("%9.6f  %9.6f  %9.6f  %9.6f\n",T[8],T[9],T[10],T[11]);
    printf("%9.6f  %9.6f  %9.6f  %9.6f\n",T[12],T[13],T[14],T[15]);
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  if(ipimfile[0]=='\0')
  {
    printf("Please specify the input image using the \"-i <image>\" flag.\n");
    exit(0);
  }

  ipim =(short *)read_nifti_image(ipimfile, &ipim_hdr);

  if(ipim  == NULL)
  {
    printf("Error reading %s, aborting ...\n", ipimfile);
    exit(0);
  }

  set_dim(ipim_dim, ipim_hdr);

  if(ipim_hdr.datatype != DT_SIGNED_SHORT && ipim_hdr.datatype != 512)
  {
    printf("\nSorry, this program only handles images of datatype\n"
    "DT_SIGNED_SHORT=4 or DT_UINT16=512. %s has datatype %d. Aborting ...\n\n", 
    ipimfile, ipim_hdr.datatype);
    free(ipim );
    exit(0);
  }

  if(opt_v)
  {
    printf("Input image = %s\n",ipimfile);
    printf("Matrix size = %d x %d x %d (voxels)\n", ipim_dim.nx, ipim_dim.ny, ipim_dim.nz);
    printf("Voxel size = %8.6f x %8.6f x %8.6f (mm3)\n", ipim_dim.dx, ipim_dim.dy, ipim_dim.dz);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  // find TPIL
  /////////////////////////////////////////////////////////////////////////////////////////////
  char orient[4]="";
  float TPIL[16];

  if(opt_v) printf("PIL transformation ...\n");
  new_PIL_transform(ipimfile,lmfile,orient, TPIL,1);
  /////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  // read PILbraincloud.nii from the $ARTHOME directory
  /////////////////////////////////////////////////////////////////////////////////////////////
  snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/PILbrain.nii",ARTHOME);

  // freed before atra() returns
  PILbraincloud = (short  *)read_nifti_image(temporaryFilename, &PILbraincloud_hdr);

  set_dim(PILbraincloud_dim, PILbraincloud_hdr);

  if(PILbraincloud==NULL)
  {
    printf("Error reading %s, aborting ...\n", temporaryFilename);
    exit(1);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////

  short *braincloud;
  int q=0;  // number of voxels in the brain mask
  braincloud = resliceImage(PILbraincloud,PILbraincloud_dim, ipim_dim,TPIL,LIN);
  for(int i=0; i<ipim_dim.nv; i++)
  {
    if(braincloud[i]<50) braincloud[i]=0; else { braincloud[i]=1; q++; }
  }
  //save_nifti_image("tt.nii", braincloud, &ipim_hdr);

  float *dispim=NULL;
  float *disp=NULL;
  int *indx;

  dispim = pixeldisp(ipim_dim, T);

  disp=(float *)calloc(q,sizeof(float));
  indx=(int *)calloc(q,sizeof(int));
  q=0;
  for(int i=0; i<ipim_dim.nv; i++)
  if(braincloud[i]==1)
  {
    disp[q++]=dispim[i]; 
  }

  hpsort(q, disp, indx);

  if(opt_v)
  {
    printf("Number of pixels inside brain mask = %d\n", q);
  }

  printf("Pixels displacement stats:\n");
  printf("min = %f\n", disp[0]);
  printf("max = %f\n", disp[q-1]);
  printf("median = %f\n", disp[q/2]);
  printf("mean = %lf\n",sample_mean(disp,q));

  free(disp);
  free(dispim);
  free(PILbraincloud);
  free(braincloud);
}

