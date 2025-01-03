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
#include <stats.h>

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
        {"-o", 1, 'o'},
        {"-h", 0, 'h'},
        {"-help", 0, 'h'},
        {0, 0, 0}
};

int opt_o=NO;

void print_help_and_exit()
{
   printf("\nUsage: one_sample_t [-h -help -b -s] -o <output image name> <input images ...>\n"
   "-b: Binarizes (0 or 1) input images before averaging\n"
   "-h or -help: Prints help message\n" 
   "-o <output image name>: Specifies the filename for the outputted average image\n\n"
   "-s: Save the output image as type short\n\n"
   "\n\n");

   exit(0);
}

float *avg(int N, char **imagefile)
{
   nifti_1_header hdr;

   int nx, ny, nz;
   float dx, dy, dz;
   int nv;
   int type;
   char *image;
   float *avg_image;

   if(N == 0) return(NULL);

   image = read_nifti_image(imagefile[0], &hdr);
   if(image==NULL) return(NULL);
   nx=hdr.dim[1]; ny=hdr.dim[2]; nz=hdr.dim[3];
   dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];

   nv = nx*ny*nz;
   avg_image = (float *)calloc(nv, sizeof(float));

   type=hdr.datatype;
   switch(type) 
   {
		case 2:
			for(int i=0; i<nv; i++) avg_image[i] = ((unsigned char *)image)[i];
			break;
		case 4:
			for(int i=0; i<nv; i++) avg_image[i] = ((short *)image)[i];
			break;
		case 8:
			for(int i=0; i<nv; i++) avg_image[i] = ((int *)image)[i];
			break;
		case 16:
			for(int i=0; i<nv; i++) avg_image[i] = ((float *)image)[i];
			break;
		case 64:
			for(int i=0; i<nv; i++) avg_image[i] = ((double *)image)[i];
			break;
		case 512:
			for(int i=0; i<nv; i++) avg_image[i] = ((unsigned short *)image)[i];
			break;
		case 256:
			for(int i=0; i<nv; i++) avg_image[i] = ((char *)image)[i];
			break;
   }

   delete image;

   for(int i=1; i<N; i++)
   {
      image = read_nifti_image(imagefile[i], &hdr);
      if(image==NULL) return(NULL);

      type=hdr.datatype; // this statement was missing; major bug fixed Aug. 8, 2012
      switch(type) {
         case 2:
            for(int j=0; j<nv; j++) avg_image[j] += ((unsigned char *)image)[j];
            break;
         case 4:
            for(int j=0; j<nv; j++) avg_image[j] += ((short *)image)[j];
            break;
         case 8:
            for(int j=0; j<nv; j++) avg_image[j] += ((int *)image)[j];
            break;
         case 16:
            for(int j=0; j<nv; j++) avg_image[j] += ((float *)image)[j];
            break;
          case 64:
            for(int j=0; j<nv; j++) avg_image[j] += ((double *)image)[j];
            break;
         case 512:
            for(int j=0; j<nv; j++) avg_image[j] += ((unsigned short *)image)[j];
            break;
         case 256:
            for(int j=0; j<nv; j++) avg_image[j] += ((char *)image)[j];
            break;
      }

      delete image;
   }

   for(int i=0; i<nv; i++) avg_image[i] /= N;

   return(avg_image);
}

float *avg4d(char *imagefile, int n)
{
   nifti_1_header hdr;

   int nx, ny, nz, nt;
   float dx, dy, dz;
   int nv;
   char *image;
   float *avg_image;

   image = read_nifti_image(imagefile, &hdr);
   if(image==NULL) 
   {
      return(NULL);
   }

   nx=hdr.dim[1]; ny=hdr.dim[2]; nz=hdr.dim[3]; nt=hdr.dim[4];
   dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];

   if(n<=0) n=nt;

   printf("Matrix size = %d x %d x %d x %d\n", nx, ny, nz, nt);
   printf("Voxel size = %f x %f x %f\n", dx, dy, dz);
   printf("Averaging %d images ...\n",n);

   if(nt == 0)
   {
      printf("nt cannot be equal to 0, aborting ...\n");
      exit(0);
   }

   nv = nx*ny*nz;
   avg_image = (float *)calloc(nv, sizeof(float));
   if(avg_image == NULL)
   {
      memory_allocation_error("\"avg_image\" in avg4d()");
   }

   for(int v=0; v<nv; v++) avg_image[v]=0.0;

   for(int i=0; i<n; i++)
   {
      switch(hdr.datatype) {
         case 2:
            for(int v=0; v<nv; v++) avg_image[v] += ((unsigned char *)image)[i*nv + v];
            break;
         case 4:
            for(int v=0; v<nv; v++) avg_image[v] += ((short *)image)[i*nv + v];
            break;
         case 8:
            for(int v=0; v<nv; v++) avg_image[v] += ((int *)image)[i*nv + v];
            break;
         case 16:
            for(int v=0; v<nv; v++) avg_image[v] += ((float *)image)[i*nv + v];
            break;
         case 64:
            for(int v=0; v<nv; v++) avg_image[v] += ((double *)image)[i*nv + v];
            break;
         case 512:
            for(int v=0; v<nv; v++) avg_image[v] += ((unsigned short *)image)[i*nv + v];
            break;
         case 256:
            for(int v=0; v<nv; v++) avg_image[v] += ((char *)image)[i*nv + v];
            break;
      }
   }

   delete image;

   for(int v=0; v<nv; v++) avg_image[v] /= n;

   return(avg_image);
}

// returns 1 if all images have the same dimensions nx, ny, and nz, 0 otherwise
int checkDimension_avgImage(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
   nifti_1_header hdr;
   short dataType;

   if(N==0) return(1);

   printf("Image %d: %s\n",1,imagefile[0]);

   hdr = read_NIFTI_hdr(imagefile[0]);
   *nx = hdr.dim[1];
   *ny = hdr.dim[2];
   *nz = hdr.dim[3];
   *dx = hdr.pixdim[1];
   *dy = hdr.pixdim[2];
   *dz = hdr.pixdim[3];
   dataType = hdr.datatype;

   for(int i=1; i<N; i++)
   {
      printf("Image %d: %s\n",i+1,imagefile[i]);
      hdr = read_NIFTI_hdr(imagefile[i]);

      if( *nx != hdr.dim[1] ||  *ny != hdr.dim[2] ||  *nz != hdr.dim[3]) 
      {
            return(0);
      }
   }

   return(1);
}

int main(int argc, char **argv)
{
   short *dfmap;
   char **imagefile;
   int n; // number of input images

   nifti_1_header hdr;
   int nx, ny, nz, nv;
   float dx, dy, dz;
   float *tmap;
   float *im;
   char outputfile[1024];
   char prefix[512];

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) 
      {
         case 'o':
            sprintf(prefix,"%s",optarg);
            opt_o=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case '?':
            print_help_and_exit();
      }
   }

   if(argc==1) print_help_and_exit();

   if(!opt_o)
   {
      printf("Please specify an output file using the -o argument.\n");
      exit(0);
   }

   n = argc-optind;

   printf("Number of input images = %d\n",n);

   if(n==0) exit(0);

   if( !checkDimension_avgImage(n, argv+optind, &nx, &ny, &nz, &dx, &dy, &dz) )
   {
      printf("\n\nAll input images must be the same size. Aborting ...");
      printf("\n\n");
      exit(0);
   }

   imagefile = (argv+optind);

   nv = nx*ny*nz;

   tmap = (float *)calloc(nv,sizeof(float));
   dfmap = (short *)calloc(nv,sizeof(short));
   im = (float *)calloc(n*nv,sizeof(float));

//hdr = read_NIFTI_hdr(imagefile[0]);
//if( hdr.dim[0] == 4 && n== 1)
//{
// write code to handle 4D images later
//}

   char *tmpim;
   for(int i=0; i<n; i++) 
   {
      tmpim = read_nifti_image(imagefile[i], &hdr);

      switch(hdr.datatype) 
      {
         case 2:
            for(int v=0; v<nv; v++) im[i*nv+v] = ((unsigned char *)tmpim)[v];
            break;
         case 4:
            for(int v=0; v<nv; v++) im[i*nv+v] = ((short *)tmpim)[v];
            break;
         case 8:
            for(int v=0; v<nv; v++) im[i*nv+v] = ((int *)tmpim)[v];
            break;
         case 16:
            for(int v=0; v<nv; v++) im[i*nv+v] = ((float *)tmpim)[v];
            break;
         case 64:
            for(int v=0; v<nv; v++) im[i*nv+v] = ((double *)tmpim)[v];
            break;
         case 512:
            for(int v=0; v<nv; v++) im[i*nv+v] = ((unsigned short *)tmpim)[v];
            break;
         case 256:
            for(int v=0; v<nv; v++) im[i*nv+v] = ((char *)tmpim)[v];
            break;
      }

      free(tmpim);
   }

   float *x;
   int df;
   double mean;
   x = (float *)calloc(n,sizeof(float));
   for(int v=0; v<nv; v++) 
   {
      for(int i=0; i<n; i++) x[i]=im[i*nv+v];
      tmap[v] = (float) one_sample_t(x, n, df, mean);
      dfmap[v]=(short)df;
   }

   hdr.datatype=16;
   sprintf(outputfile,"%s_tmap.nii",prefix);
   save_nifti_image(outputfile, tmap, &hdr);

   hdr.datatype=4;
   sprintf(outputfile,"%s_df.nii",prefix);
   save_nifti_image(outputfile, dfmap, &hdr);

   free(tmap);
   free(dfmap);
   free(im);
   free(x);

   return 0;
}
