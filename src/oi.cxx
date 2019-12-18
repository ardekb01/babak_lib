#include <babak_lib.h>
#include <sph.h>
#include <minmax.h>
#include <math.h>

/////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-o",1,'o'},
   {"-i",1,'i'},
   {"-n",1,'n'},
   {"-v",0,'v'},
   {0, 0,  0}
};

int main(int argc, char **argv)
{
   float *im;
   int n=2; // default
   char imfile[512]="";
   nifti_1_header hdr;

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'n':
            n = atoi(optarg);
            if(n<2) n=2;
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'i':
            sprintf(imfile,"%s",optarg);
            break;
         case '?':
            exit(0);
      }
   }

//////////////////////////////////////////////////////////////////////////////////////
   int nx, ny, nz, np, nv;
   float dx, dy, dz;
   float min=0.0, max=0.0;
   float avgsize=0.0;
   double oi=0.0;

   if( imfile[0]=='\0' )
   {
      printf("Please specify an input image using -i option.\n");
      exit(0);
   }

   im = (float *)read_nifti_image(imfile, &hdr);
   nx = hdr.dim[1];
   ny = hdr.dim[2];
   nz = hdr.dim[3];
   dx = hdr.pixdim[1];
   dy = hdr.pixdim[2];
   dz = hdr.pixdim[3];
   nv = nx*ny*nz;
   np = nx*ny;

   minmax(im,nv,min,max);

   if(opt_v)
   {
      printf("Input Image = %s\n",imfile);
      printf("Matrix size = %d x %d x %d\n", nx, ny, nz);
      printf("Voxel size = %f x %f x %f\n", dx, dy, dz);
      printf("Min=%f, Max=%f\n",min,max);
      printf("n=%d\n", n);
   }

   for(int i=0; i<nv; i++)
   {
      im[i] /= max;
      avgsize += im[i];
   }

   minmax(im,nv,min,max);

   // Note: ln(x)/ln(2) = log_2(x)
   // In any case the ln(2)'s cancel so we can just use ln() instead of log_2()
   for(int i=0; i<nv; i++)
   if( im[i]>0.0 )
   {
      // oi += im[i] * log( (double) im[i])/log(2.0);
      oi += im[i] * log( (double) im[i]);
   }

   // oi = 1 + log(2.0)*oi/(avgsize*log(1.0*n));
   oi = 1 + oi/(avgsize*log(1.0*n));

   printf("Average size=%7.2f, OI=%7.5lf\n", avgsize, oi);
}
