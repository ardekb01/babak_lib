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
#include <spm_analyze.h>
#include <babak_lib.h>
#include <f2c.h>

#define YES 1
#define NO 0

int opt;

static struct CmdOption options[] =
{
  {"-version", 0, 'V'},
  {"-V", 0, 'V'},
  {"-verbose", 0, 'v'},
  {"-v", 0, 'v'},
  {"-float", 0, 'f'},
  {"-f", 0, 'f'},
  {"-short", 0, 's'},
  {"-s", 0, 's'},
  {"-prefix", 1, 'p'},
  {"-p", 1, 'p'},
  {"-scale", 1, 'c'},
  {0, 0, 0}
};

int opt_f=NO;
int opt_s=NO;
int opt_p=NO;
int opt_scale=NO;
	
void print_help_and_exit()
{
  printf("\nUsage: scaleImage [-V -v -f -s] [-scale <scale>] [-prefix <prefix>] <image files ...>\n\n"
  "-V or -version : prints date of last revision\n\n"
  "-v or -verbose : sets verbose mode\n\n"
  "-f or -float : output image is saved as type 'float'\n\n"
  "-s or -short : output image is saved as type 'short'\n\n"
  "-p or -prefix : prefix for naming output images (default='scl_')\n\n"
  "-scale : each voxel value is multiplied by this number (default=1.0)\n\n"
  );
  exit(0);
}

int main(int argc, char **argv)
{
  nifti_1_header hdr;
  char prefix[512];
  char tempfile[512];
  char opfile[1024];
  char *im;
  char **imagefile;
  float scale=1.0;
  float dx,dy,dz;
  float *fim;
  int n;
  int nx, ny, nz, datatype;
  int nv;
  short *sim;
  unsigned char *ucim;

  while( (opt=getoption(argc, argv, options)) != -1)
  {
    switch (opt) {
      case 'V':
        printf("\n\nLast modified: April 15, 2005\n\n");
        exit(0);
      case 'v':
        opt_v=YES;
        break;
      case 's':
        opt_s=YES;
        break;
      case 'f':
        opt_f=YES;
        break;
      case 'p':
        snprintf(prefix,sizeof(prefix),"%s",optArg);
        opt_p=YES;
        break;
      case 'c':
        scale = atof(optArg);
        opt_scale=YES;
        break;
      case '?':
        print_help_and_exit();
    }
  }

  if(argc==1) print_help_and_exit();

  if(opt_s && opt_f) 
  {
    printf("\n -s (-short) and -f (-float) options cannot be used together!\n");
    exit(0);
  }

  if(!opt_p) snprintf(prefix,sizeof(prefix),"scl_");

  if(opt_v && opt_p) printf("\nOutput image prefix is '%s'\n",prefix);

  if(opt_v) printf("\nScale = %f\n",scale);

  n = argc-optind;
  if(opt_v) printf("\nNumber of input images = %d\n",n);

  if(n<=0) exit(0);

  imagefile = argv+optind;

  if(opt_v) 
  {
    printf("\nInput images:\n");
    for(int i=0; i<n; i++)
      printf("\tImage %d = %s\n",i+1,imagefile[i]);
  }

  for(int i=0; i<n; i++)
  {
    im = NULL;
    im = read_nifti_image(imagefile[i], &hdr);
    nx = hdr.dim[1]; ny = hdr.dim[2]; nz = hdr.dim[3];
    dx = hdr.pixdim[1]; dy = hdr.pixdim[2]; dz = hdr.pixdim[3];
    datatype = hdr.datatype;

    nv = nx*ny*nz;

    // I added this because some SPM ANALYZE images have negative dimensions
    if(dx<0) dx *= -1; 
    if(dy<0) dy *= -1; 
    if(dz<0) dz *= -1;

    if(im==NULL) continue;

    // if(datatype!=4 && datatype!=16)
    if(datatype!=2 && datatype!=4 && datatype!=16)
    {
      printf("\nCannot process datatype %d, skipping this image ...\n",datatype);
      continue;
    }

    if(opt_v)
    {
      printf("\nProcessing image %d\n", i+1);
      printf("\tImage filename = %s\n",imagefile[i]);
      if( datatype == 2 ) printf("\tInput image data type = 'unsigned char'\n");
      if( datatype == 4 ) printf("\tInput image data type = 'short'\n");
      if( datatype == 16 ) printf("\tInput image data type = 'float'\n");
      printf("\tMatrix size = %d x %d x %d\n",nx,ny,nz);
      printf("\tVoxel size = %f x %f x %f\n",dx,dy,dz);
    }

    getfilename(tempfile,imagefile[i]);
    snprintf(opfile,sizeof(opfile),"%s%s",prefix,tempfile);

    if(opt_v) printf("\tOutput image filename = %s\n",opfile);

    if( datatype==2 )
    {
      ucim = (unsigned char *)im;
      for(int v=0; v<nv; v++)  ucim[v] = (unsigned char)( ucim[v]*scale );
      if(opt_v) printf("\tOutput image data type = 'unsigned char'\n");
      save_nifti_image(opfile, ucim, &hdr);
    }

    if( datatype==4 )
    {
      sim = (short *)im;

      if(opt_f)
      {
        fim = (float *)calloc(nv,sizeof(float));
        for(int v=0; v<nv; v++)  fim[v] = sim[v]*scale;
        if(opt_v) printf("\tOutput image data type = 'float'\n");
        save_nifti_image(opfile, fim, &hdr);
        free(fim);
      }
      else 
      {
        for(int v=0; v<nv; v++)  sim[v] = (short)( sim[v]*scale + 0.5);
        if(opt_v) printf("\tOutput image data type = 'short'\n");
        save_nifti_image(opfile, sim, &hdr);
      }
    }

    if( datatype==16 )
    {
      fim = (float *)im;

      if(opt_s)
      {
        sim = (short *)calloc(nv,sizeof(short));
        for(int v=0; v<nv; v++)  sim[v] = (short) (fim[v]*scale + 0.5);
        if(opt_v) printf("\tOutput image data type = 'short'\n");
        save_nifti_image(opfile, sim, &hdr);
        free(sim);
      }
      else
      {
        for(int v=0; v<nv; v++)  fim[v] *= scale;
        if(opt_v) printf("\tOutput image data type = 'float'\n");
        save_nifti_image(opfile, fim, &hdr);
      }
    }

    free(im);
  }
}
