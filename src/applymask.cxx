#include <stdlib.h>
#include <malloc.h>  
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <unistd.h>
//#include "spm_analyze.h"
#include "babak_lib.h"

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
  {"-m", 1, 'm'},
  {"-i", 1, 'i'},
  {"-o", 1, 'o'},
  {"-v", 0, 'v'},
  {0, 0, 0}
};

void print_help_and_exit()
{
  printf("\nUsage: applymask -m <mask> -i <input image> -o <output image>\n\n");
  exit(0);
}

short *readmask(char *maskfile,int &nx_mask,int &ny_mask,int &nz_mask)
{
  nifti_1_header hdr;
  float dx, dy, dz;
  int nv;
  int type;
  char *image;
  short *mask;

  if(opt_v)
  {
    printf("Reading mask file %s ...\n",maskfile);
  }

  image = read_nifti_image(maskfile, &hdr);

  if(image==NULL) {
    printf("Reading maskfile %s failed, aborting ...\n",maskfile);
    exit(0);
  }

  nx_mask = hdr.dim[1]; ny_mask = hdr.dim[2]; nz_mask = hdr.dim[3];
  dx = hdr.pixdim[1]; dy = hdr.pixdim[2]; dz = hdr.pixdim[3];
  type = hdr.datatype;

  if(opt_v)
  {
    printf("\tMatrix size: %d x %d x %d\n",nx_mask, ny_mask, nz_mask);
    printf("\tVoxel size: %6.4f x %6.4f x %6.4f\n",dx, dy, dz);
    printf("\tData type: %d\n",type);
  }

  nv = nx_mask*ny_mask*nz_mask;

  mask = (short *)calloc(nv, sizeof(short));

  if(type==512) type=4;
  if(type==256) type=2;
  if(type==768) type=8;

  switch(type) {
    case 2:
      for(int i=0; i<nv; i++) mask[i] = ((unsigned char *)image)[i];
      break;
    case 4:
      for(int i=0; i<nv; i++) mask[i] = ((short *)image)[i];
      break;
    case 8:
      for(int i=0; i<nv; i++) mask[i] = ((int *)image)[i];
      break;
    case 16:
      for(int i=0; i<nv; i++) mask[i] = (short)((float *)image)[i];
      break;
  }

  delete image;

  return(mask);
}

int main(int argc, char **argv)
{
  nifti_1_header hdr;
  short *mask;
  int nx_mask, ny_mask, nz_mask;
  int nx, ny, nz;
  int nv;
  float dx, dy, dz;
  char outputfile[1024]="";
  char inputfile[1024]="";
  char maskfile[1024]="";
  int type;
  short *image;

  while( (opt=getoption(argc, argv, options)) != -1)
  {
    switch (opt) {
    case 'm':
      sprintf(maskfile,"%s",optarg);
      break;
    case 'i':
      sprintf(inputfile,"%s",optarg);
      break;
    case 'o':
      sprintf(outputfile,"%s",optarg);
      break;
    case 'v':
      opt_v=YES;
      break;
    case '?':
      print_help_and_exit();
    }
  }

  if(argc==1) print_help_and_exit();

  if( maskfile[0]=='\0' )
  {
    printf("Please specify a mask file using the -m option. Aborting ...\n");
    exit(0);
  }

  if(opt_v)
  {
    printf("Mask file = %s\n",maskfile);
  }

  if( inputfile[0]=='\0' )
  {
    printf("Please specify an input image using the -i option. Aborting ...\n");
    exit(0);
  }

  if(opt_v)
  {
    printf("Input image = %s\n",inputfile);
  }

  if( outputfile[0]=='\0' )
  {
    printf("Please specify an output image using the -o option. Aborting ...\n");
    exit(0);
  }

  if(opt_v)
  {
    printf("Output image = %s\n",outputfile);
  }

  nx_mask = ny_mask = nz_mask = 0;
  mask=readmask(maskfile,nx_mask,ny_mask,nz_mask);
  nv = nx_mask * ny_mask * nz_mask;

  if(opt_v)
  {
    printf("Reading input image %s ...\n",inputfile);
  }
  image = (short *)read_nifti_image(inputfile, &hdr);

  if(image==NULL) {
    printf("Reading input image %s failed, aborting ...\n",inputfile);
    exit(0);
  }

  nx = hdr.dim[1]; ny = hdr.dim[2]; nz = hdr.dim[3];
  dx = hdr.pixdim[1]; dy = hdr.pixdim[2]; dz = hdr.pixdim[3];
  type = hdr.datatype;

  if(opt_v)
  {
    printf("\tMatrix size: %d x %d x %d\n",nx, ny, nz);
    printf("\tVoxel size: %6.4f x %6.4f x %6.4f\n",dx, dy, dz);
    printf("\tData type: %d\n",type);
  }

  if(nx!= nx_mask || ny!=ny_mask || nz!=nz_mask)
  {
    printf("Mask dimensions do not match input image dimensions, aborting ...\n");
  }

  for(int v=0; v<nv; v++) if(mask[v]==0) image[v]=0;

  if(opt_v)
  {
    printf("Writing output image %s ...\n",outputfile);
  }
  save_nifti_image(outputfile, image, &hdr);

  delete image;
  delete mask;
}

