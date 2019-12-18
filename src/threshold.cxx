#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "spm_analyze.h"
#include "babak_lib.h"
#include <f2c.h>

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
  {"-v", 0, 'v'},
  {"-cc", 1, 'c'},
  {"-t1", 1, '1'},
  {"-t2", 1, '2'},
  {"-o", 1, 'o'},
  {"-i", 1, 'i'},
  {0, 0, 0}
};

void print_help_and_exit()
{
	exit(0);
}

char opt_t1=NO;
char opt_t2=NO;

int main(int argc, char **argv)
{
  char *im;
  int ccthresh=0;
  float thresh1=0, thresh2=0;
  char op_image_file[1024]="";
  int n1=0; 
  int n2=0;
  float *ip_image;
  float *op_image;

  nifti_1_header hdr;

	char inputfile[512];

	int nx, ny, nz, nv;
	float dx, dy, dz;

  while( (opt=getoption(argc, argv, options)) != -1)
  {
    switch (opt) {
      case 'v':
        opt_v=YES;
        break;
      case 'c':
        ccthresh=atoi(optarg);
        break;
      case 'i':
        sprintf(inputfile,"%s",optarg);
        break;
      case 'o':
        sprintf(op_image_file,"%s",optarg);
        break;
      case '1':
        opt_t1=YES;
        thresh1 = atof(optarg);
        break;
      case '2':
        opt_t2=YES;
        thresh2 = atof(optarg);
        break;
      case '?':
        print_help_and_exit();
    }
  }

  if(argc==1) print_help_and_exit();

  if(opt_v)
  {
    printf("CC threshold = %d\n",ccthresh);
    if(opt_t1) printf("Thresh 1 = %f\n",thresh1);
    if(opt_t2) printf("Thresh 2 = %f\n",thresh2);
    if(op_image_file[0] != '\0') printf("Output image = %s\n", op_image_file);
  }
   
  ip_image = (float *)read_nifti_image(inputfile, &hdr);
  nx = hdr.dim[1]; ny = hdr.dim[2]; nz = hdr.dim[3];
  dx = hdr.pixdim[1]; dy = hdr.pixdim[2]; dz = hdr.pixdim[3];

  nv = nx*ny*nz;
  op_image = (float *)calloc(nv, sizeof(float));
  im = (char *)calloc(nv, sizeof(char));

  for(int i=0; i<nv; i++) im[i]=0;

  n1=0;
  if(opt_t1)
  {
    for(int i=0; i<nv; i++)
    {
      if(ip_image[i]>thresh1)
      {
        im[i]=1;
        n1++;
      }
    }
  }

  n2=0;
  if(opt_t2)
  {
    for(int i=0; i<nv; i++)
    {
      if(ip_image[i]<thresh2)
      {
        im[i]=1;
        n2++;
      }
    }
  }

  if(ccthresh>0) 
  {
    int ncc, ntcc;

    thr_Connected_Component(im, ccthresh, nx, ny, nz, &ncc, &ntcc);
    printf("\n%d clusters of %d had size >= %d\n",ntcc,ncc,ccthresh);
  }

  for(int i=0; i<nv; i++) if(im[i]==1) op_image[i]=ip_image[i];

  if(op_image_file[0] != '\0') save_nifti_image(op_image_file, op_image, &hdr);

  if(opt_t1) printf("%d\n",n1);
  if(opt_t2) printf("%d\n",n2);

  free(op_image);
  free(im);
}
